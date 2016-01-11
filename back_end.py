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
from os import remove, stat
from shutil import move
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


    def NCBI_search(self):
        """
        Submit search to NCBI and return the records.
        """
        handle = Entrez.esearch(db=self.database, term=self.term,
                                usehistory="y", retmax=100000000)
        record = Entrez.read(handle)
        handle.close()

        return record


    def Record_processor(self, record):
        """
        Splits the record returned by Entrez into sparate variables and returns
        them.
        """
        count = int(record["Count"]) # Int
        IDs = record["IdList"] # List
        webenv = record["WebEnv"] # String
        query_key = record["QueryKey"] #String

        assert count == len(IDs)

        return count, IDs, webenv, query_key


    def NCBI_history_fetch(self, count, IDs, webenv, query_key, Bsize, Run):
        """
        Fetches results from NCBI using history.
        """
        try:
            a = open(self.outfile, 'r')
            a.close()
        except:
            a = open(self.outfile, 'w')
            a.close()
        if Run == 1 and stat(self.outfile).st_size != 0:
            self.ReDownloader(IDs, webenv, query_key, Bsize)
        else:

            if count == 0 and self.gui == 0:
                sys.exit("No records found in database!")
            elif count == 0:
                self.no_match.emit("No sequences in the database matched your " "query.")
                return None

            else:
                outfile = open(self.outfile, 'a')
                if self.gui == 1:
                    self.max_seq.emit(count)
                if Bsize > count:
                    Bsize = count
                for start in range(0, count, Bsize):
                    if start + Bsize < count:
                        end = start + Bsize
                    else:
                        end = count
                    print("Downloading record %i to %i of %i" %(start+1, end,
                                                                count))

                    if self.gui == 1:
                        self.prog_data.emit(end)

                    # Make sure that the program carries on despite server "hammering" errors.
                    attempt = 0
                    while True:
                        try:
                            fetch_handle = Entrez.efetch(db=self.database,
                                                         rettype="fasta",
                                                         retstart=start,
                                                         retmax=Bsize,
                                                         webenv=webenv,
                                                         query_key=query_key)
                            data = fetch_handle.read()
                            fetch_handle.close()
                            if data.startswith("<?"):
                                raise ValueError("NCBI server error.")
                            else:
                                data = data.replace("\n\n","\n")
                                break
                        except ValueError:
                            if attempt < 5:
                                print("NCBI is retuning XML instead of sequence"
                                      " data. Trying the same chunk again in "
                                      "8\'\'.")
                                attempt += 1
                                sleep(8)
                                pass
                            else:
                                print("Too many errors in a row. Let's make a "
                                      "larger 20\'\' pause and try again.")
                                attempt = 0
                                sleep(20)
                                pass
                    outfile.write(data)
        try:
            outfile.close()
        except:
            pass
        self.ReDownloader(IDs, webenv, query_key, Bsize)


    def ReDownloader(self, IDs, webenv, query_key, Bsize):
        """
        Checks for missing sequences.
        """
        print("Checking for sequences that did not download... Please wait.")
        ver_IDs = self.Error_finder(self.outfile)
        missing_IDs = []
        for i in IDs:
            if i not in ver_IDs:
                missing_IDs.append(i)
        numb_missing = len(missing_IDs)
        IDs = missing_IDs # Improve performance on subsequent runs
        if numb_missing == 0:
            print("All sequences were downloaded correctly. Good!")
            if self.gui == 0:
                sys.exit("Program finished without error.")
        else:
            print("%s sequences did not download correctly (or at all). "
                  "Retrying..." %(numb_missing))
            self.NCBI_history_fetch(numb_missing, IDs, webenv, query_key, Bsize,
                                    2)


    def Error_finder(self, target_file):
        """
        Looks for errors in the output fasta and retruns a list of necessary retries.
        """
        target_handle = open(target_file, 'r')
        verified_IDs = set()

        for lines in target_handle:
            if lines.startswith(">"):
                ID = re.search("gi\|.*?\|", lines).group(0)[3:-1]
                verified_IDs.add(ID)

        target_handle.close()
        return verified_IDs


    def runEverything(self):
        """
        Run the functions in order.
        """
        batch_size = 2000
        Entrez.email = self.email

        rec = self.NCBI_search()
        count, IDs, webenv, query_key = self.Record_processor(rec)
        self.NCBI_history_fetch(count, IDs, webenv, query_key, batch_size, 1)


def main():
    """
    Main function. Defines how the arguments get passed to the rest of the
    program.
    """
    dl = Downloader(sys.argv[1], sys.argv[3], sys.argv[2], sys.argv[4], 0)
    dl.runEverything()


if __name__ == '__main__':
    main()
