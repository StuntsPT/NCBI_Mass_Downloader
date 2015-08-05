#NCBI Mass Sequence Downloader

##Description:
This program will download sequences *en masse* from several NCBI databases (at the user's chioce).
It will run on both python2 and python3!

##Installation
Not required. Just clone the repository or uncompress one of the releases and run from there.

##Usage:
To use the GUI version:

    python NCBI_downloader.py

To use the command line version:

    python NCBI_downloader.py "user@email-address.com" "Query term" "database" outfile.fasta

Notes:
* The program will not overwrite anything in the output file, but rather append sequences to it;
* For now the *fasta* format is hard-coded in the program, but this may change should there be demand for it;
* *Query term* can take any argument just like in the website search engine (eg. "Lacerta monticola[organism]");
* The program uses the Entrez module from [biopython](https://github.com/biopython/biopython). A big Kudos to the authors.

##Dependencies:
* pyhton (2 or 3) (and the standard lib);
* [PyQt5](http://www.riverbankcomputing.com/software/pyqt/intro) for the GUI version.

##License:
GPLv3

##Known limitations:
* There is **NO** argument checking of any kind. If there is a problem with your arguments, you will see no *hand holding* of any kind, just plain python error messages.

##Found a bug?
Or would like a feature added? Or maybe drop some feedback?
Just open an issue on github.
