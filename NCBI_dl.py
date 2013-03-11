#!/usb/bin/python2

#Usage: NCBI_dl.py "Query term" "database" outfile.fasta
from Bio import Entrez
from sys import argv

Entrez.email = "f.pinamartins@gmail.com"

handle = Entrez.esearch(db=argv[2],term=argv[1],retmax=1000000)
record = Entrez.read(handle)
IDs = ",".join(record["IdList"])


handle2 = Entrez.efetch(db=argv[2], id=IDs, rettype="fasta")

outfile = open(argv[3],'w')

outfile.write(handle2.read())

outfile.close()
