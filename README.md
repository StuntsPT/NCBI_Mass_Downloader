#NCBI Mass Sequence Downloader

##Description:
This program will download sequences *en masse* from several NCBI databases (at the user's chioce).

##Usage:
To use the GUI version:

    python3 NCBI_downloader.py

To use the command line version:

    python3 NCBI_downloader.py "user@email-address.com" "Query term" "database" outfile.fasta

Notes:
* The program will not overwrite anything in the output file, but rather append sequences to it!
* For now the *fasta* format is hard-coded in the program, but this may change at a later date.
* *Query term* can take any argument just like in the website search engine (eg. "Lacerta monticola[organism]").

##Dependencies:
* pyhton2 (and the standard lib);
* [biopython](https://github.com/biopython/biopython);
* [PyQt4](http://www.riverbankcomputing.com/software/pyqt/intro) for the GUI version;

##License:
GPLv3

##Known limitations:
There is **NO** argument checking of any kind. If there is a problem with your arguments, you will see no *hand holding* of any kind.
Just plain error messages.

##Found a bug?
Or would like a feature added? Or maybe just wanto to drop some feedback?
Just open an issue on github!