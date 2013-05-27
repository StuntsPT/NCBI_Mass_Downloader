#NCBI Mass Sequence Downloader

##Description:
This program will download sequences "*en masse*" from several NCBI databases (at the user's chioce).

##Usage:
    python2 NCBI_dl.py "user@email-address.com" "Query term" "database" outfile.fasta
Note that the program will not overwrite anything in the output file, but rather append sequences to it!
For now the *fasta* format is hard-coded in the program, but this may change at a later date.

Just use the program with all the arguments and let the program download your sequences.
*Query term* can take any argument just like in the website search engine (eg. "Lacerta monticola[organism]").

##Dependencies:
* pyhton2 (and the standard lib);
* [biopython](https://github.com/biopython/biopython);

##License:
GPLv3

##Found a bug?
Or would like a feature added? Or maybe just wanto to drop some feedback?
Just open an issue on github!