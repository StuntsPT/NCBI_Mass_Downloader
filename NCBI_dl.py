#!/usb/bin/python2

#Usage: NCBI_dl.py "Query term" "database" outfile.fasta
from Bio import Entrez
from sys import argv

Entrez.email = "f.pinamartins@gmail.com"

handle = Entrez.esearch(db=argv[2],term=argv[1],retmax=1000000)
record = Entrez.read(handle)
IDs = record["IdList"]


fasta_data = []

for i in IDs:
	handle2 = Entrez.efetch(db=argv[2], id=i, rettype="fasta")
	fasta_data.append(handle2.read())

outfile = open(argv[3],'w')

for i in fasta_data:
	outfile.write(str(i))

outfile.close()
