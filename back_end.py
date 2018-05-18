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
        super(Downloader, self).__init__()
        self.run = 0
        self.webenv = ""
        self.query_key = ""

    def ncbi_search(self, database, term):
        """
        Submit search to NCBI and return the records.
        """
        handle = Entrez.esearch(db=database, term=term, usehistory="y")
        record = Entrez.read(handle)
        handle.close()

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

    def main_organizer(self, count, ids, b_size):
        """
        Defines what tasks need to be performed, handles NCBI server errors and
        writes the sequences into the outfile.
        """
        try:
            if self.run == 1 and stat(self.outfile).st_size != 0:
                self.re_downloader(ids, b_size, count)
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

                if ids == []:
                    fetch_func = self.fetch_by_history
                    fetch_args = start, b_size
                else:
                    fetch_func = self.fetch_by_id
                    fetch_args = ids, b_size
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
            self.re_downloader(ids, b_size, count)

    def re_downloader(self, ids, b_size, count):
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
            if len(ver_ids) == count:
                print("All sequences were downloaded correctly. Good!")
                if self.gui == 0:
                    sys.exit("Program finished without error.")
                else:
                    self.finished.emit("Download finished successfully!")
            else:
                if ids == []:
                    missing_ids = self.id_fetcher(count, done_set=ver_ids)
                else:
                    missing_ids = []
                    for i in ids:
                        if i not in ver_ids:
                            missing_ids.append(i)
                numb_missing = len(missing_ids)
                ids = missing_ids  # Improve performance on subsequent runs
                print("%s sequences did not download correctly (or at all). "
                      "Retrying..." % (numb_missing))
                self.run = 2
                self.main_organizer(numb_missing, ids, b_size)

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

    def fetch_by_id(self, ids, b_size):
        """
        Fetches NCBI data based on the IDs, rather than a search query. Returns
        the data handle string.
        """
        id_handle = Entrez.efetch(db=self.database,
                                  id=ids,
                                  rettype="fasta",
                                  retmode="text",
                                  retmax=b_size)
        data = id_handle.read()
        id_handle.close()

        return data

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

    def id_fetcher(self, count, done_set=None):
        """
        Gets a list of all ACC numbers matching a query and returns them as a
        list.
        """
        chunk = 10000  # 10000 is the maximum retmax allowed for ACC idtype
        ids = []

        for i in range(0, count, chunk):
            success = False
            if i > count:
                i = count
            while not success:
                try:
                    iter_handle = Entrez.efetch(db=self.database,
                                                webenv=self.webenv,
                                                query_key=self.query_key,
                                                retmax=chunk,
                                                rettype="acc", retstart=i)
                    success = True
                except HTTPError:
                    print("Got an HTTPError. Let's wait 8'' and try again.")
                    sleep(8)
            if done_set is None:
                ids += [x.rstrip() for x in iter_handle]
            else:
                ids += [x.rstrip() for x in iter_handle if x not in done_set]
            iter_handle.close()

        return ids

    def run_everything(self):
        """
        Run the functions in order.
        """
        batch_size = 3000
        Entrez.email = "frmartins@ciencias.ulisboa.pt"
        Entrez.api_key = "bbceccfdf97b6b7e06e93c918e010f1ecf09"
        self.run = 1

        rec = self.ncbi_search(self.database, self.term)
        try:
            count, self.webenv, self.query_key = self.record_processor(rec)
            ids = []
        except TypeError:
            return None
        if self.database == "genome":
            ids = self.translate_genome(count)
            count = len(ids)
            self.database = "nucleotide"
            self.run = 2
        self.main_organizer(count, ids, batch_size)


def main():
    """
    Main function. Defines how the arguments get passed to the rest of the
    program.
    """
    dler = Downloader(sys.argv[2], sys.argv[1], sys.argv[3], 0)
    dler.run_everything()


if __name__ == '__main__':
    main()
