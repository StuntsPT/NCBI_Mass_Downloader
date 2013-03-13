#!/usb/bin/python2

#Usage: NCBI_dl.py "Query term" "database" outfile.fasta
from Bio import Entrez
from sys import argv

#Set global vars:
user_email = "f.pinamartins@gmail.com"
database = argv[2]
search_term = argv[1]
output_file = argv[3]
batch_size = 1000


def NCBI_search(Uemail, DB, ST):
    #Submit search to NCBI and return the records
    Entrez.email = Uemail
    handle = Entrez.esearch(db=DB,term=ST,usehistory="y",retmax=10000000)
    record = Entrez.read(handle)
    handle.close()

    return record

def Record_processor(record):
    #Processes the record into sparate usefull information
    count = int(record["Count"])
    IDs = record["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]

    assert count == len(IDs)

    return count, IDs, webenv, query_key

def NCBI_fetch(output_file, count, IDs, webenv, query_key, Bsize):
    #Fetches results from NCBI
    outfile = open(output_file,'w')
    if Bsize > count:
        Bsize = count
    for start in range(0,count,Bsize):
        end = start + Bsize
        print("Downloading record %i to %i of %i") % (start+1, end, count)
        fetch_handle = Entrez.efetch(db=database, rettype="fasta", retstart=start, retmax=Bsize, webenv=webenv, query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        outfile.write(data)

    outfile.close()

#Run everything
rec = NCBI_search(user_email, database, search_term)
count, IDs, webenv, query_key = Record_processor(rec)
NCBI_fetch(output_file, count, IDs, webenv, query_key, batch_size)
