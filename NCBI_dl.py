#!/usb/bin/python2

#Usage: NCBI_dl.py "Query term" "database" outfile.fasta
from Bio import Entrez
from sys import argv

Entrez.email = "f.pinamartins@gmail.com"

handle = Entrez.esearch(db=argv[2],term=argv[1],usehistory="y",retmax=1000000)
record = Entrez.read(handle)
handle.close()

count = int(record["Count"])
IDs = record["IdList"]
webenv = record["WebEnv"]
query_key = record["QueryKey"]

print(count)
print(len(IDs))

assert count == len(IDs)

batch_size = 1000

outfile = open(argv[3],'w')
for start in range(0,count,batch_size):
    end = min(count,batch_size)
    print("Going to download record %i to %i") % (start+1, end)
    fetch_handle = Entrez.efetch(db=argv[2], rettype="fasta", restart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
    data = fetch_handle.read()
    fetch_handle.close()
    outfile.write(data)

outfile.close()
