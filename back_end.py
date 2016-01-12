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

# Usage: back_end.py "user@email-address.com" "Query term" "database" outfile.fasta

import sys
import re
from os import stat
from time import sleep

import Entrez

class Downloader(object):
    def __init__(self, email, database, term, outfile, gui):
        self.email = email
        self.database = database
        self.term = term
        self.outfile = outfile
        self.gui = gui
        super(Downloader, self).__init__()


    def ncbi_search(self, database, term):
        """
        Submit search to NCBI and return the records.
        """
        handle = Entrez.esearch(db=database, term=term, usehistory="y",
                                retmax=100000000)
        record = Entrez.read(handle)
        handle.close()

        return record


    def record_processor(self, record):
        """
        Splits the record returned by Entrez into sparate variables and returns
        them.
        """
        count = int(record["Count"]) # Int
        IDs = record["IdList"] # List
        webenv = record["WebEnv"] # String
        query_key = record["QueryKey"] #String

        assert count == len(IDs)

        if count == 0 and self.gui == 0:
            sys.exit("Your serch query returned no results!")

        elif count == 0:
            self.no_match.emit("Your serch query returned no results!")
            return None

        return count, IDs, webenv, query_key


    def main_organizer(self, count, IDs, webenv, query_key, b_size, Run):
        """
        Defines what tasks need to be performed, handles NCBI server errors and
        writes the sequences into the outfile.
        """
        try:
            if Run == 1 and stat(self.outfile).st_size != 0:
                gui_end = self.re_downloader(IDs, webenv, query_key, b_size)
                return None
        except OSError:
            pass

        outfile = open(self.outfile, 'a')
        if b_size > count:
            b_size = count
        for start in range(0, count, b_size):
            if start + b_size < count:
                end = start + b_size
            else:
                end = count
            print("Downloading record %i to %i of %i" %(start+1, end, count))

            if self.gui == 1:
                self.max_seq.emit(count)
                self.prog_data.emit(end)

            if Run == 1:
                fetch_func = self.fetch_by_history
                fetch_args = start, b_size, webenv, query_key
            else:
                fetch_func = self.fetch_by_id
                fetch_args = IDs, b_size
            # Make sure that the program carries on despite server "hammering" errors.
            attempt = 0
            while True:
                try:
                    data = fetch_func(*fetch_args)
                    if data.startswith("<?"):
                        raise ValueError("NCBI server error.")
                    else:
                        data = data.replace("\n\n", "\n")
                        break
                except:
                    if attempt < 5:
                        print("NCBI is retuning XML instead of sequence data. "
                              "Trying the same chunk again in 8\'\'.")
                        attempt += 1
                        sleep(8)
                    else:
                        print("Too many errors in a row. Let's make a larger "
                              "20\'\' pause and try again.")
                        attempt = 0
                        sleep(20)
            outfile.write(data)

        outfile.close()
        self.re_downloader(IDs, webenv, query_key, b_size)


    def re_downloader(self, IDs, webenv, query_key, b_size):
        """
        Checks for missing sequences.
        """
        print("Checking for sequences that did not download... Please wait.")
        ver_ids = self.error_finder(self.outfile)
        missing_ids = []
        for i in IDs:
            if i not in ver_ids:
                missing_ids.append(i)
        numb_missing = len(missing_ids)
        IDs = missing_ids # Improve performance on subsequent runs
        if numb_missing == 0:
            print("All sequences were downloaded correctly. Good!")
            if self.gui == 0:
                sys.exit("Program finished without error.")

        else:
            print("%s sequences did not download correctly (or at all). "
                  "Retrying..." %(numb_missing))
            self.main_organizer(numb_missing, IDs, webenv, query_key, b_size, 2)


    def error_finder(self, target_file):
        """
        Looks for errors in the output fasta and retruns a list of necessary retries.
        """
        target_handle = open(target_file, 'r')
        verified_ids = set()

        for lines in target_handle:
            if lines.startswith(">"):
                ID = re.search("gi\|.*?\|", lines).group(0)[3:-1]
                verified_ids.add(ID)

        target_handle.close()
        return verified_ids


    def fetch_by_id(self, IDs, b_size):
        """
        Fetches NCBI data based on the IDs, rather than a search query. Returns
        the data handle string.
        """
        id_handle = Entrez.efetch(db=self.database,
                                  id=IDs,
                                  rettype="fasta",
                                  retmode="text",
                                  retmax=b_size)
        data = id_handle.read()
        id_handle.close()

        return data


    def fetch_by_history(self, start, b_size, webenv, query_key):
        """
        Fetches NCBI data based on the provided search query. Returns the data
        handle string.
        """
        hist_handle = Entrez.efetch(db=self.database,
                                    retstart=start,
                                    rettype="fasta",
                                    retmode="text",
                                    retmax=b_size,
                                    webenv=webenv,
                                    query_key=query_key)
        data = hist_handle.read()
        hist_handle.close()

        return data


    def run_everything(self):
        """
        Run the functions in order.
        """
        batch_size = 3000
        Entrez.email = self.email

        rec = self.ncbi_search(self.database, self.term)
        try:
            count, IDs, webenv, query_key = self.record_processor(rec)
        except TypeError:
            return None
        self.main_organizer(count, IDs, webenv, query_key, batch_size, 1)


def main():
    """
    Main function. Defines how the arguments get passed to the rest of the
    program.
    """
    dler = Downloader(sys.argv[1], sys.argv[3], sys.argv[2], sys.argv[4], 0)
    dler.run_everything()


if __name__ == '__main__':
    main()
