#!/usb/bin/python2

#Usage: NCBI_dl.py "Query term" "database" outfile.fasta
from Bio import Entrez
from sys import argv
from shutil import move
from os import remove
import re

#Set global vars:
user_email = "f.pinamartins@gmail.com"
database = argv[2]
search_term = argv[1]
output_file = argv[3]
batch_size = 1000

Entrez.email = user_email

def NCBI_search(ST):
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

def NCBI_history_fetch(output_file, count, IDs, webenv, query_key, Bsize):
    #Fetches results from NCBI using history
    outfile = open(output_file,'a')
    if Bsize > count:
        Bsize = count
    for start in range(0,count,Bsize):
        if start + Bsize < count:
            end = start + Bsize
        else:
            end = count
        print("Downloading record %i to %i of %i") % (start+1, end, count)
        #TODO: Make it a "try" to enable retries on server errors
        fetch_handle = Entrez.efetch(db=database, rettype="fasta", retstart=start, retmax=Bsize, webenv=webenv, query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        outfile.write(data)

    outfile.close()
    ReDownloader(output_file, IDs)

def ReDownloader(output_file, IDs):
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
        print("%s sequences did not download correctly. Retrying...") %(len(missing_IDs))
        count, IDs, webenv, query_key = NCBI_post(IDs)
        NCBI_history_fetch(output_file, count, IDs, webenv, query_key, 1000)


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
        elif lines..strip().startswith("<") or lines.startswith("\n"):
            pass
        else:
            new_file.write(lines)

    original_file.close()
    new_file.close()
    remove(temp_file)

    return verified_IDs



#Run everything
rec = NCBI_search(search_term)
count, IDs, webenv, query_key = Record_processor(rec)
NCBI_history_fetch(output_file, count, IDs, webenv, query_key, batch_size)