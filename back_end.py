#!/usr/bin/python3
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

# Usage: back_end.py "Query term" "database" outfile.fasta

import sys
import re
import tempfile
import pickle
from os import stat
from time import sleep
from json import decoder

import requests


class ProgramDone(Exception):
    """
    Simple class to raise a custom error (terminated by user)
    """
    pass


class Downloader():
    def __init__(self, database, term, outfile, gui):
        self.database = database
        self.term = term
        self.outfile = outfile
        self.gui = gui
        super(Downloader, self).__init__()
        self.api_key = "bbceccfdf97b6b7e06e93c918e010f1ecf09"
        self.original_count = 0
        self.terminated = False
        self.accn_cache = tempfile.TemporaryFile()


    def ncbi_search(self, database, term, history="y", retmax=0, retstart=0):
        """
        Submit search to NCBI and return the WebEnv & QueryKey.
        """
        if self.terminated is True:
            return
        record = {}
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {"db": database,
                         "term": term,
                         "retmode": "json",
                         "usehistory": history,
                         "idtype": "acc",
                         "retmax": retmax,
                         "retstart": retstart,
                         "api_key": self.api_key}

        attempt = 0
        while attempt < 5:
            try:
                handle = requests.get(url, params=search_params, timeout=30)
                record["qkey"] = [handle.json()["esearchresult"]["querykey"]]
                record["webenv"] = handle.json()["esearchresult"]["webenv"]
                record["count"] = int(handle.json()["esearchresult"]["count"])
                record["accn"] = handle.json()["esearchresult"]["idlist"]
                break
            except KeyError:
                err = handle.text
                self.finish(False,
                            "NCBI returned an error:\n" + err + "\nExiting.")
            except:
                print("Got an empty reply or a Timeout from NCBI."
                      " Let's wait 8'' and try again.")
                attempt += 1
                sleep(8)
        else:
            self.finish(False, "Too many errors in a row.\n"
                        "Please wait a few minutes and try "
                        "again.")

        if record["count"] == 0 and self.gui == 0:
            sys.exit("Your serch query returned no results!")
        elif record["count"] == 0 and self.gui == 1:
            self.no_match.emit("Your serch query returned no results!")
            return None

        return record

    def use_efetch(self, start, b_size, webenv, query_key):
        """
        Fetches NCBI data based on the provided search query or acession
        numbers. Returns the fasta string.
        """
        if self.terminated is True:
            raise ProgramDone
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {"db": self.database,
                         "retmode": "text",
                         "rettype": "fasta",
                         "api_key": self.api_key,
                         "WebEnv": webenv,
                         "query_key": query_key,
                         "retstart": start,
                         "retmax": b_size,
                         "term": self.term}
        attempt = 0
        while attempt < 5:
            try:
                handle = requests.get(url, params=fetch_params, timeout=30)
                fasta = handle.text
                if not fasta.startswith(">"):
                    raise ValueError("NCBI server error.")
                fasta = fasta.replace("\n\n", "\n")
                break
            except:
                attempt += 1
                print("NCBI is not retuning sequence data or we're getting a "
                      "Timeout. Trying the same chunk again in 8\'\'.")
                sleep(8)
        else:
            self.finish(False, "Too many errors in a row.\n"
                        "Please wait a few minutes and try "
                        "again.")

        return fasta


    def use_epost(self, accession_numbers, webenv):
        """
        Uploads data to an 'artificial' history bank on NCBI's servers.
        Returns a WebEnv and a query_key
        """
        if self.terminated is True:
            raise ProgramDone
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi"
        post_params = {"db": self.database,
                       "api_key": self.api_key}
        if webenv is not None:
            post_params["WebEnv"] = webenv
        attempt = 0
        while attempt < 5:
            try:
                handle = requests.post(url,
                                       params=post_params,
                                       data={"id": accession_numbers},
                                       timeout=30)
                break
            except:
                attempt += 1
                print("NCBI is either retuning an error or we're getting "
                      "a Timeout. Trying the same chunk again in 8\'\'.")
                sleep(8)
        else:
            self.finish(False, "Too many errors in a row.\n"
                        "Please wait a few minutes and try "
                        "again.")
        if webenv is None:
            webenv = re.search("<WebEnv>.*</WebEnv>",
                               handle.text).group()[8:-9]
        query_key = re.search("<QueryKey>.*</QueryKey>",
                              handle.text).group()[10:-11]

        return webenv, query_key


    def main_organizer(self, count, b_size, query_key, webenv):
        """
        Defines what tasks need to be performed, handles NCBI server errors and
        writes the sequences into the outfile.
        """
        try:
            stat(self.outfile).st_size
            file_exists = True
        except OSError:
            missing_accns = None
            file_exists = False
        except TypeError:
            return
        if file_exists:
            if stat(self.outfile).st_size != 0:
                missing_accns = self.missing_checker()
                count = len(missing_accns)
            else:
                missing_accns = None


        if missing_accns is not None:
            self.artificial_history(missing_accns)
        else:
            self.actual_downloader(count, b_size, query_key, webenv)


    def missing_checker(self):
        """
        Checks if any sequences did not download correctlly.
        If discrepancies between the downloaded sequences and the Accesions are
        found, a new missing sequences list is generated and an
        attempt is made to fetch the missing sequences.
        """
        print("Checking for sequences that did not download... Please wait.")
        ver_ids = self.fasta_parser(self.outfile)
        if self.accn_cache.tell() != 0:
            # If accecsion numbers are cached, just load them instead of
            # downloading them again
            self.accn_cache.seek(0)
            ncbi_accn_set = pickle.load(self.accn_cache)
            print("Using cached accession numbers.")
        else:
            retmax = 50000
            if self.original_count <= retmax:
                retmax = self.original_count
            ncbi_accn_set = set()
            for i in range(0, self.original_count, retmax):
                if i + retmax < self.original_count:
                    end = i + retmax
                else:
                    end = self.original_count
                print("Downloading accession %i to %i of "
                      "%i" % (i + 1, end, self.original_count))

                subset = set(self.ncbi_search(self.database,
                                              self.term,
                                              "n",
                                              retmax=retmax,
                                              retstart=i)["accn"])
                ncbi_accn_set = ncbi_accn_set.union(subset)

            # Remove any Master records from the accn set:
            # See https://www.biostars.org/p/305310/#305317
            master_records = set()
            for accn in ncbi_accn_set:
                if bool(re.search('[A-Z]{4}0+(\.\d){0,}$', accn)):
                    master_records.add(accn)
                    print("WARNING: Master record found and "
                          "removed: %s." % (accn))
            ncbi_accn_set = ncbi_accn_set - master_records

            # Create an accecsion number cache. This should avoid subsequent
            # accession number downloads.
            pickle.dump(ncbi_accn_set, self.accn_cache, pickle.HIGHEST_PROTOCOL)


        if ver_ids == ncbi_accn_set:
            self.finish(success=True)

        missing_ids = ncbi_accn_set - ver_ids

        return missing_ids


    def fasta_parser(self, target_file):
        """
        Parses a FASTA file and retruns a set of found Accession numbers
        """
        target_handle = open(target_file, 'r')
        verified_ids = set()

        for lines in target_handle:
            if lines.startswith(">"):
                seqid = re.match("([^\s]+)", lines).group(0)[1:]
                verified_ids.add(seqid)

        target_handle.close()
        return verified_ids


    def artificial_history(self, accns):
        """
        Takes a list of accn numbers, posts it to NCBI to generate an
        'artificial' history (webenv and query_key). This avoids download by
        accn, and works around the protein issue described here:
        https://www.biostars.org/p/476638/
        """
        if self.terminated is True:
            raise ProgramDone
        outfile = open(self.outfile, 'a')
        max_batch = 200
        if len(accns) < max_batch:
            max_batch = len(accns)
        webenv = None
        accn_list = list(accns)
        print("Creating an aritificial 'history' to work around direct "
              "accession number download.")
        for i in range(0, len(accns), max_batch):
            if self.terminated is True:
                return
            if i + max_batch < len(accns):
                end = i + max_batch
            else:
                end = len(accns)

            accns_str = ",".join(accn_list[i:end])

            # Use epost
            print("Uploading accessions %i to %i of %i" % (i + 1,
                                                           end,
                                                           len(accns)))
            webenv, query_key = self.use_epost(accns_str, webenv)

            # Use efetch!
            print("Downloading sequences %i to %i of %i" % (i + 1,
                                                            end,
                                                            len(accns)))
            if self.gui == 1:
                self.max_seq.emit(len(accns))
                self.prog_data.emit(end)

            fasta_data = self.use_efetch(0,
                                         max_batch,
                                         webenv=webenv,
                                         query_key=query_key)
            outfile.write(fasta_data)

        outfile.close()
        self.missing_checker()


    def genome_deprecation(self):
        """
        Emits a message on how to obtian genomes from NCBI.
        """
        message = ("NCBI has deprecated the 'Genomes' dataset from e-utils. "
                  "The cannonical way to get genome data is currently to use "
                  "their 'Datasets' utility:\n"
                  "https://www.ncbi.nlm.nih.gov/datasets/genomes/\n\n"
                  "For more information please visit 'NCBI Insights':\n"
                  "https://ncbiinsights.ncbi.nlm.nih.gov/2020/09/10/genomic-data/")

        self.finish(False, message)
        return 0


    def finish(self, success, msg=""):
        """
        Run this funcion whenever downloads are finished.
        """
        if success is True:
            print("All sequences were downloaded correctly. Good!")
            if self.gui == 0:
                sys.exit("Program finished without error.")
            else:
                self.finished.emit("Download finished successfully!")
                raise ProgramDone

        else:
            if self.gui == 0:
                sys.exit("Program finished with some failures.\n" + msg)
            else:
                self.finished.emit("Program finished with some "
                                   "failures.\n" + msg)
                raise ProgramDone


    def actual_downloader(self, count, b_size, query_key, webenv):
        """
        Manages downloads
        """
        outfile = open(self.outfile, 'a')
        if b_size > count:
            b_size = count
        for start in range(0, count, b_size):
            if start + b_size < count:
                end = start + b_size
            else:
                end = count
            print("Downloading record %i to %i of %i" % (start + 1, end, count))

            if self.gui == 1:
                self.max_seq.emit(count)
                self.prog_data.emit(end)

            data = self.use_efetch(start,
                                   b_size,
                                   webenv=webenv,
                                   query_key=query_key)
            outfile.write(data)

        outfile.close()
        self.main_organizer(count, b_size, query_key, webenv)


    def run_everything(self):
        """
        Run the functions in order.
        """
        try:
            if self.database == "genome":
                self.genome_deprecation()
                return

            batch_size = 200

            record = self.ncbi_search(self.database, self.term)
            count = record["count"]
            self.original_count = count

            self.main_organizer(count, batch_size, record["qkey"], record["webenv"])
        except ProgramDone:
            return

def main():
    """
    Main function. Defines how the arguments get passed to the rest of the
    program.
    """
    dler = Downloader(sys.argv[2], sys.argv[1], sys.argv[3], 0)
    dler.run_everything()


if __name__ == '__main__':
    main()
