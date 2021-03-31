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
import requests
import json
from os import stat
from time import sleep

import Entrez


class Downloader(object):
    def __init__(self, database, term, outfile, gui):
        self.database = database
        self.term = term
        self.outfile = outfile
        self.gui = gui
        self.terminated = False
        self.failures = [[], 0]
        self.retry_threshold = 5
        super(Downloader, self).__init__()
        self.run = 0
        self.webenv = ""
        self.query_key = ""

    def ncbi_search(self, database, term, history="y", retmax=0):
        """
        Submit search to NCBI and return the WebEnv & QueryKey.
        """
        record = {}
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {"db": database,
                         "term": term,
                         "retmode": "json",
                         "usehistory": history,
                         "idtype": "acc",
                         "retmax": retmax}

        handle = requests.get(url, params=search_params)
        record["qkey"] = handle.json()["esearchresult"]["querykey"]
        record["webenv"] = handle.json()["esearchresult"]["webenv"]
        record["count"] = int(handle.json()["esearchresult"]["count"])
        record["accn"] = handle.json()["esearchresult"]["idlist"]

        return record


    def record_processor(self, record):
        """
        Splits the record returned by Entrez into sparate variables and returns
        them.
        """
        count = int(record["Count"])  # Int
        webenv = record["WebEnv"]  # String
        query_key = record["QueryKey"]  # String

        if count == 0 and self.gui == 0:
            sys.exit("Your serch query returned no results!")

        elif count == 0:
            self.no_match.emit("Your serch query returned no results!")
            return None

        return count, webenv, query_key

    def main_organizer(self, count, b_size):
        """
        Defines what tasks need to be performed, handles NCBI server errors and
        writes the sequences into the outfile.
        """
        try:
            if self.run == 1 and stat(self.outfile).st_size != 0:
                self.re_downloader(b_size, count)
                return None
        except OSError:
            pass

        outfile = open(self.outfile, 'a')
        if b_size > count:
            b_size = count
        for start in range(0, count, b_size):
            if self.terminated is False:
                if start + b_size < count:
                    end = start + b_size
                else:
                    end = count
                print("Downloading record %i to %i of %i" % (start + 1, end,
                                                             count))

                if self.gui == 1:
                    self.max_seq.emit(count)
                    self.prog_data.emit(end)


                fetch_func = self.fetch_by_history
                fetch_args = start, b_size

                # Make sure that the program carries on despite server
                # "hammering" errors.
                attempt = 0
                while self.terminated is False:
                    try:
                        data = fetch_func(*fetch_args)
                        if data.startswith("<?"):
                            raise ValueError("NCBI server error.")
                        else:
                            data = data.replace("\n\n", "\n")
                            break
                    except ValueError:
                        if attempt < 5:
                            print("NCBI is retuning XML instead of sequence "
                                  "data. Trying the same chunk again in 8\'\'.")
                            attempt += 1
                            sleep(8)
                        else:
                            print("Too many errors in a row. Let's make a "
                                  " larger 20\'\' pause and try again.")
                            attempt = 0
                            sleep(20)
                if self.terminated is False:
                    outfile.write(data)

        outfile.close()
        if self.terminated is False:
            self.re_downloader(b_size, count)

    def re_downloader(self, b_size, count):
        """
        Checks if any sequences did not download correctlly.
        If discrepancies between the downloaded sequences and the IdList are
        found, a new IdList is generated, based on what is missing and an
        attempt is made to fetch the missing sequences.
        """
        if self.terminated is True:
            return
        else:
            print("Checking for sequences that did not download... Please "
                  "wait.")
            ver_ids = self.fasta_parser(self.outfile)
            print("Getting an accession number list to compare with the local "
                  "file... Please wait")
            ncbi_accn_set = set(self.ncbi_search(self.database,
                                                 self.term, "n", count)["accn"])
            # Remove any Master records from the accn set:
            master_records = set()
            for accn in ncbi_accn_set:
                if bool(re.search('[A-Z]{4}0+(\.\d){0,}$', accn)):
                    master_records.add(accn)
                    print("WARNING: Master record found and "
                          "removed: %s." % (accn))
            ncbi_accn_set = ncbi_accn_set - master_records

            if ver_ids == ncbi_accn_set:
                print("All sequences were downloaded correctly. Good!")
                if self.gui == 0:
                    sys.exit("Program finished without error.")
                else:
                    self.finished.emit("Download finished successfully!")
            else:
                missing_ids = ncbi_accn_set - ver_ids

                if self.failures[0] != missing_ids:
                    self.failures[0] = missing_ids
                    self.failures[1] = 0
                else:
                    self.failures[1] += 1
                numb_missing = len(missing_ids)
                if numb_missing == 0:
                    print("All sequences were downloaded correctly. Good!")
                    if self.gui == 0:
                        sys.exit("Program finished without error.")
                    else:
                        self.finished.emit("Download finished successfully!")

                elif self.failures[1] < self.retry_threshold:
                    print("%s sequences did not download correctly (or at all)."
                          " Retrying..." % (numb_missing))
                    self.run = 2
                    self.main_organizer(numb_missing, b_size)
                else:
                    print("NOTICE: After %s retries, not all sequences were "
                          "downloaded correctly.=-(" % (self.retry_threshold))
                    print("A list of failed downloads can be found on %s.failed"
                          % (self.outfile))
                    fail_log = open(self.outfile + ".failed", "w")
                    for i in missing_ids:
                        fail_log.write(i + "\n")
                    fail_log.close()
                    if self.gui == 0:
                        sys.exit("Program finished without error.")
                    else:
                        self.finished.emit("Download finished with some "
                                           "failures!\nPlease check the file "
                                           "%s.failed for a detailed "
                                           "list." % (self.outfile))


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

    def fetch_by_accn(self, accns, b_size):
        """
        Fetches NCBI data based on the IDs, rather than a search query. Returns
        the data handle string.
        """
        # id_handle = Entrez.efetch(db=self.database,
        #                           id=ids,
        #                           rettype="fasta",
        #                           retmode="text",
        #                           retmax=b_size)
        # data = id_handle.read()
        # id_handle.close()

        return

    def fetch_by_history(self, start, b_size):
        """
        Fetches NCBI data based on the provided search query. Returns the data
        handle string.
        """
        hist_handle = Entrez.efetch(db=self.database,
                                    retstart=start,
                                    rettype="fasta",
                                    retmode="text",
                                    retmax=b_size,
                                    webenv=self.webenv,
                                    query_key=self.query_key)
        data = hist_handle.read()
        hist_handle.close()

        return data

    def translate_genome(self, count):
        """
        Translates genome query IDs into a nucleotide query IDs, since NCBI has
        deprecated the use of the "genome" database, and the old genome IDs.
        http://www.ncbi.nlm.nih.gov/books/NBK25499/
        """
        import urllib

        ids = self.id_fetcher(count)

        nuc_acc_list = []
        query_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + \
                    "elink.fcgi?dbfrom=genome&db=nucleotide&id="
        for genome_id in ids:
            tmplist = []
            xml = urllib.request.urlopen(query_url + genome_id)
            for content in xml:
                if content.endswith(b"</Id>\n"):
                    gid = re.search("<Id>.*</Id>",
                                    content.decode('utf-8')).group()[4:-5]
                    tmplist.append(gid)
            nuc_acc_list += tmplist[1:]

        return nuc_acc_list


    def run_everything(self):
        """
        Run the functions in order.
        """
        batch_size = 3000
        Entrez.email = "frmartins@ciencias.ulisboa.pt"
        Entrez.api_key = "bbceccfdf97b6b7e06e93c918e010f1ecf09"
        self.run = 1

        record = self.ncbi_search(self.database, self.term)
        count = record["count"]
        self.query_key = record["qkey"]
        self.webenv = record["webenv"]

        if self.database == "genome":
            ids = self.translate_genome(count)
            count = len(ids)
            self.database = "nucleotide"
            self.run = 2
        self.main_organizer(count, batch_size)


def main():
    """
    Main function. Defines how the arguments get passed to the rest of the
    program.
    """
    dler = Downloader(sys.argv[2], sys.argv[1], sys.argv[3], 0)
    dler.run_everything()


if __name__ == '__main__':
    main()
