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
        Processes the record into sparate usefull information
        """
        count = int(record["Count"])
        IDs = record["IdList"]
        webenv = record["WebEnv"]
        query_key = record["QueryKey"]

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
            self.ReDownloader(IDs, webenv, query_key)
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

                    # Make sure that even on server errors the program carries on.
                    # If the servers are dead, well, you were not going anywhere anyway...
                    while True:
                        try:
                            fetch_handle = Entrez.efetch(db=self.database,
                                                         rettype="fasta",
                                                         retstart=start,
                                                         retmax=Bsize,
                                                         webenv=webenv,
                                                         query_key=query_key)
                            break
                        except:
                            pass
                    data = fetch_handle.read()
                    fetch_handle.close()
                    outfile.write(data)
        try:
            outfile.close()
        except:
            pass
        self.ReDownloader(IDs, webenv, query_key)


    def ReDownloader(self, IDs, webenv, query_key):
        """
        Checks for missing sequences.
        """
        print("Checking for sequences that did not download... Please wait.")
        ver_IDs = self.Error_finder()
        missing_IDs = set()
        for i in IDs:
            if i not in ver_IDs:
                missing_IDs.add(i)
        IDs = missing_IDs # Improve performance on subsequent runs
        if len(missing_IDs) == 0:
            print("All sequences were downloaded correctly. Good!")
            if self.gui == 0:
                sys.exit("Program finished without error.")
        else:
            print("%s sequences did not download correctly (or at all). "
                  "Retrying..." %(len(missing_IDs)))
            self.NCBI_history_fetch(len(IDs), IDs, webenv, query_key, 1000, 2)


    def Error_finder(self):
        """
        Looks for errors in the output fasta and retruns a list of necessary retries.
        """
        temp_file = self.outfile + ".tmp"
        move(self.outfile, temp_file)
        original_file = open(temp_file, 'r')
        new_file = open(self.outfile, 'w')
        verified_IDs = set()

        for lines in original_file:
            # Why do this? Sometimes, empty lines are downloaded, and this
            # method removes them. It also removes some "^<" that ocasinally
            # show up. Not sure whay that is, but this is a good (albeit slow)
            # workaround.
            if lines.startswith(">"):
                ID = re.search("gi\|.*?\|", lines).group(0)[3:-1]
                verified_IDs.add(ID)
                new_file.write(lines)
            elif lines.strip().startswith("<") or lines.startswith("\n"):
                pass
            else:
                new_file.write(lines)

        original_file.close()
        new_file.close()
        remove(temp_file)

        return verified_IDs


    def runEverything(self):
        """
        Run the functions in order.
        """
        batch_size = 1000
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
