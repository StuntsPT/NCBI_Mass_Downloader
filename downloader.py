#!/usb/bin/python3
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

#Usage: NCBI_dl.py "user@email-address.com" "Query term" "database" outfile.fasta
#Note that the program will not overwrite the outfile, but rather append sequences to it!
#Meant to be used from NCBI_dl.py, but can also be used "stand-alone".
from Bio import Entrez
from sys import argv
from shutil import move
from os import remove, stat
import re

def NCBI_search(ST, database):
    #Submit search to NCBI and return the records
    handle = Entrez.esearch(db=database,term=ST,usehistory="y",retmax=10000000)
    record = Entrez.read(handle)
    handle.close()

    return record

def NCBI_post(IDs):
    #Submit id_list to NCBI via epost and return the records
    IDs_string = ",".join(IDs)
    handle = Entrez.epost(database,id=IDs_string,retmax=10000000)
    record = Entrez.read(handle)
    handle.close()

    count = len(IDs)
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]

    return count, IDs, webenv, query_key

def Record_processor(record):
    #Processes the record into sparate usefull information
    count = int(record["Count"])
    IDs = record["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]

    assert count == len(IDs)

    return count, IDs, webenv, query_key

def NCBI_history_fetch(output_file, count, IDs, webenv, query_key, Bsize, Run, database):
    #Fetches results from NCBI using history
    try:
        a = open(output_file,'r')
        a.close()
    except:
        a = open(output_file,'w')
        a.close()
    if Run == 1 and stat(output_file).st_size != 0:
        ReDownloader(output_file, IDs, database)
    else:
        outfile = open(output_file,'a')
        if Bsize > count:
            Bsize = count
        for start in range(0,count,Bsize):
            if start + Bsize < count:
                end = start + Bsize
            else:
                end = count
            print("Downloading record %i to %i of %i") % (start+1, end, count)
            #Make sure that even on server errors the program carries on.
            #If the servers are dead, well, you were not going anywhere anyway...
            while True:
                try:
                    fetch_handle = Entrez.efetch(db=database, rettype="fasta", retstart=start, retmax=Bsize, webenv=webenv, query_key=query_key)
                    break
                except:
                    print(1)
                    pass
            data = fetch_handle.read()
            fetch_handle.close()
            outfile.write(data)

    outfile.close()
    ReDownloader(output_file, IDs, database)

def ReDownloader(output_file, IDs, database):
    #Check for missing sequences:
    print("Checking for sequences that did not download... Please wait.")
    ver_IDs = Error_finder(output_file)
    missing_IDs = set()
    for i in IDs:
        if i not in ver_IDs:
            missing_IDs.add(i)
    IDs = missing_IDs #Improve performance on subsequent runs
    if len(missing_IDs) == 0:
        print("All sequences were downloaded correctly. Good!")
        quit("Program finished without error.")
    else:
        print("%s sequences did not download correctly (or at all). Retrying...") %(len(missing_IDs))
        count, IDs, webenv, query_key = NCBI_post(IDs)
        NCBI_history_fetch(output_file, count, IDs, webenv, query_key, 1000, 2, database)


def Error_finder(output_file):
    #Looks for errors in the output fasta and retruns a list of necessary retries
    temp_file = output_file + ".tmp"
    move(output_file, temp_file)
    original_file = open(temp_file,'r')
    new_file = open(output_file,'w')
    verified_IDs = set()

    for lines in original_file:
        if lines.startswith(">"):
            ID = re.search("gi\|.*?\|",lines).group(0)[3:-1]
            verified_IDs.add(ID)
            new_file.write("\n" + lines) #TODO: remove first empty line from file
        elif lines.strip().startswith("<") or lines.startswith("\n"):
            pass
        else:
            new_file.write(lines)

    original_file.close()
    new_file.close()
    remove(temp_file)

    return verified_IDs

def main():
    #Set vars or get them from the GUI:
    user_email = argv[1]
    database = argv[3]
    search_term = argv[2]
    output_file = argv[4]

    runEverything(user_email, database, search_term, output_file)

def runEverything(user_email, database, search_term, output_file):
    #Run the functions in order
    batch_size = 1000
    Entrez.email = user_email

    rec = NCBI_search(search_term, database)
    count, IDs, webenv, query_key = Record_processor(rec)
    NCBI_history_fetch(output_file, count, IDs, webenv, query_key, batch_size, 1, database)

if __name__ == '__main__':
    main()