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
from multiprocessing.dummy import Pool
from itertools import repeat
from random import randrange
from os import stat
from time import sleep
from json import decoder

#from icecream import ic

import requests


class ProgramDone(Exception):
    """
    Simple class to raise a custom error (terminated by user)
    """
    pass


class Downloader():
    def __init__(self, database, term, outfile, verification, gui):
        self.database = database
        self.term = term
        self.outfile = outfile
        self.gui = gui
        super(Downloader, self).__init__()
        self.api_key = "bbceccfdf97b6b7e06e93c918e010f1ecf09"
        self.original_count = 0
        self.terminated = False
        self.accn_cache = tempfile.TemporaryFile()
        self.batch_size = 200
        self.progress = 0
        self.pool_size = 8  # Max. 10 threads allowed by NCBI using an API key)
        self.verification = verification
        self.verification_attempt = 0


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

        sleep(1)
        return record


    def use_efetch(self, start, webenv, query_key, count=0, epost_batch=0):
        """
        Fetches NCBI data based on the provided search query or acession
        numbers. Returns the fasta string.
        """
        if self.terminated is True:
            raise ProgramDone

        # Ensure workers start async
        sleep(0.01 * randrange(0, 250, 5))

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {"db": self.database,
                         "retmode": "text",
                         "rettype": "fasta",
                         "api_key": self.api_key,
                         "WebEnv": webenv,
                         "query_key": query_key,
                         "retstart": start,
                         "retmax": self.batch_size,
                         "term": self.term}


        if epost_batch != 0:
            print_start = epost_batch + 1
            end = epost_batch + self.batch_size
        else:
            print_start = start + 1
            end = start + self.batch_size
        if end > count:
            end = count

        print("Downloading record %i to %i of %i" % (print_start, end, count))
        if self.gui == 1:
            if end > self.progress:
                self.progress = end
                self.prog_data.emit(end)

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

        print("Uploading %i accession numbers to NCBI's history "
              "server." % (self.batch_size))

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


    def main_organizer(self, count, query_key, webenv, missing_accns=None):
        """
        Defines what tasks need to be performed, handles NCBI server errors and
        writes the sequences into the outfile.
        """
        if self.verification is True:
            try:
                stat(self.outfile).st_size
                file_exists = True
            except OSError:
                file_exists = False
            except TypeError:
                return
            if file_exists:
                if stat(self.outfile).st_size != 0 and missing_accns is None:
                    missing_accns = self.missing_checker()
                    count = len(missing_accns)

            if missing_accns is not None:
                missing_accns = self.artificial_history(missing_accns)
            else:
                self.actual_downloader(count, query_key, webenv)
                missing_accns = self.missing_checker()
            count = len(missing_accns)
            self.main_organizer(count, query_key, webenv, missing_accns)
        else:
            with open(self.outfile, 'w') as otf:
                otf.close()
            self.actual_downloader(count, query_key, webenv)
            self.finish(success=True)


    def missing_checker(self):
        """
        Checks if any sequences did not download correctlly.
        If discrepancies between the downloaded sequences and the Accesions are
        found, a new missing sequences list is generated and an
        attempt is made to fetch the missing sequences.
        """
        self.verification_attempt += 1
        if self.verification_attempt >= 6:
            self.finish(False, "After 5 failed attempts to verify the download,"
                        " it is apparent that some accession numbers cannot be "
                        "matched to the FASTA titles. If you can perform manual"
                        " validation, consider using the '-nv' switch to skip "
                        " the automatic verification step.")
        print("Checking for sequences that did not download... Please wait.")
        ver_ids = self.fasta_parser(self.outfile)
        if self.accn_cache.tell() != 0:
            # If accecsion numbers are cached, just load them instead of
            # downloading them again
            self.accn_cache.seek(0)
            ncbi_accn_set = pickle.load(self.accn_cache)
            print("Using cached accession numbers.")
        else:
            # with open('data.pickle', 'rb') as f:
            # # The protocol version used is detected automatically, so we do not
            # # have to specify it.
            #     ncbi_accn_set = pickle.load(f)
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
            # with open('data.pickle', 'wb') as f:
            #     # Pickle the 'data' dictionary using the highest protocol available.
            #     pickle.dump(ncbi_accn_set, f, pickle.HIGHEST_PROTOCOL)

        missing_ids = ncbi_accn_set - ver_ids

        if missing_ids != set():
            not_missing = self.check_unconformant(missing_ids, ver_ids)
            missing_ids = missing_ids - not_missing

        # debug_file = open("missing_ids.txt", "w")
        # for line in missing_ids:
        #     debug_file.write(line + "\n")
        # debug_file.close()

        if missing_ids == set():
            self.finish(success=True)

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


    def check_unconformant(self, not_found, local_set):
        """
        Makes sure that IDs that were not found in the local FASTA are not
        simply encoded in a non-standard way:
        eg "SOMETHING|Accession|SOMETHING"
        Returns a set of missing IDs with any matched entries removed
        """
        not_missing = set()
        for title in local_set:
            if "|" in title:
                not_missing.add(re.search("\|.*\|", title).group()[1:-1])
                not_missing.add(re.search("\|.*$", title).group()[1:].replace("|", ""))
                # not_missing.add(re.search("\|.*$", title).group()[1:].replace("|", "") + "+")
                # not_missing.add(re.sub(".$", r"+\g<0>", re.search("\|.*$", title).group()[1:].replace("|", "")))
                # not_missing.add(re.search("\|.*$", title).group()[1:].replace("|", "_"))

        not_missing = not_missing.intersection(not_found)

        # between_pipes = {re.search("\|.*\|", x).group()[1:-1]
        #                  if "|" in x else "" for x in local_set}
        # cut_by_pipe = {re.search("\|.* ", x).group()[1:-1].replace("|", "")
        #                if "|" in x else "" for x in local_set}
        # replaced_pipe = {re.search("\|.* ", x).group()[1:-1].replace("|", "_")
        #                  if "|" in x else "" for x in local_set}
        #
        # not_missing_bp = between_pipes.intersection(not_found)
        # not_missing_cbp = cut_by_pipe.intersection(not_found)
        # not_missing_rp = replaced_pipe.intersection(not_found)
        #
        # not_missing = not_missing_bp.add(not_missing_cbp).add(not_missing_rp)

        return not_missing


    def artificial_history(self, accns):
        """
        Takes a list of accn numbers, posts it to NCBI to generate an
        'artificial' history (webenv and query_key). This avoids download by
        accn, and works around the protein issue described here:
        https://www.biostars.org/p/476638/
        """
        if self.terminated is True:
            raise ProgramDone

        count = len(accns)

        # Split the accn list into a list of 200 accns strings
        accn_strs = self.splitter(list(accns), self.batch_size, "s")
        batches = [200 * x for x in range(0, len(accn_strs))]

        accn_containers = self.splitter(accn_strs, 500)
        batch_containers = self.splitter(batches, 500)

        if self.gui == 1:
            self.max_seq.emit(len(accns))
            self.prog_data.emit(0)

        print("Creating an aritificial 'history'...")
        # Every 100K sequences, we force the outfile closed, and start a new
        # worker pool, in order to force periodic filesystem writes.
        for batch, accns_slice in zip(batch_containers, accn_containers):
            pool = Pool(self.pool_size)
            outfile = open(self.outfile, 'a')
            for fasta in pool.starmap(self.generate_and_get_from_history,
                                      zip(accns_slice,
                                          repeat(None),
                                          repeat(count),
                                          batch)):

                outfile.write(fasta)

            outfile.close()

        missing = self.missing_checker()

        return missing


    def generate_and_get_from_history(self, accns_str, webenv, count, batch):
        """
        Generates the artificail history, and immediately gets the results
        from it.
        """
        webenv, query_key = self.use_epost(accns_str, webenv)
        fasta_data = self.use_efetch(0, webenv, query_key, count, batch)

        return fasta_data


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


    def actual_downloader(self, count, query_key, webenv):
        """
        Manages downloads
        """
        if self.gui == 1:
            self.max_seq.emit(count)
            self.prog_data.emit(0)

        if count % self.batch_size != 0:
            batches = list(range(0, count, self.batch_size))
        else:
            batches = list(range(0, count, self.batch_size))[:-1]

        # Split the batches list into a list of lists, each with 100K sequences
        # (500 batches each)
        containers = self.splitter(batches, 500)

        # Every 100K sequences, we force the outfile closed, and start a new
        # worker pool, in order to force periodic filesystem writes.
        for batch in containers:
            pool = Pool(self.pool_size)
            outfile = open(self.outfile, 'a')
            for fasta in pool.starmap(self.use_efetch,
                                      zip(batch,
                                          repeat(webenv),
                                          repeat(query_key),
                                          repeat(count))):

                outfile.write(fasta)

            outfile.close()

        return


    def splitter(self, lts, size, res="l"):
        """
        Splits a "list to split" (lts) into a list of lists with length 'size'
        and returns it. If 'res' == 's', return a list of strings instead
        (with list elements joined by ',')
        Based on https://stackoverflow.com/a/2215676/3091595
        """
        if res == "l":
            new_list = [lts[i:i + size] for i in range(0, len(lts), size)]
        elif res == "s":
            new_list = [",".join(lts[i:i + size])
                        for i in range(0, len(lts), size)]

        return new_list


    def run_everything(self):
        """
        Run the functions in order.
        """
        try:
            if self.database == "genome":
                self.genome_deprecation()
                return

            record = self.ncbi_search(self.database, self.term)
            count = record["count"]
            self.original_count = count

            self.main_organizer(count, record["qkey"], record["webenv"])
        except ProgramDone:
            return

def main():
    """
    Main function. Defines how the arguments get passed to the rest of the
    program.
    """
    dler = Downloader(sys.argv[2], sys.argv[1], sys.argv[3], sys.argv[4], 0)
    dler.run_everything()


if __name__ == '__main__':
    main()
