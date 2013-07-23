#NCBI Mass Sequence Downloader

##Description:
This program will download sequences *en masse* from several NCBI databases (at the user's chioce).
It will run on both python2 and python3!

##Usage:
To use the GUI version:

    python NCBI_downloader.py

To use the command line version:

    python NCBI_downloader.py "user@email-address.com" "Query term" "database" outfile.fasta

Notes:
* The program will not overwrite anything in the output file, but rather append sequences to it!
* For now the *fasta* format is hard-coded in the program, but this may change at a later date.
* *Query term* can take any argument just like in the website search engine (eg. "Lacerta monticola[organism]").
* The program uses the Entrez module from [biopython](https://github.com/biopython/biopython). A big Kudos to the authors;

##Dependencies:
* pyhton (2 or 3) (and the standard lib);
* [PyQt4](http://www.riverbankcomputing.com/software/pyqt/intro) for the GUI version;

##License:
GPLv3

##Known limitations:
* There is **NO** argument checking of any kind. If there is a problem with your arguments, you will see no *hand holding* of any kind.
Just plain error messages.
* There is a bug in Biopython 1.61 with urllib when being used with python 3. This bug was corrected with [this](https://github.com/biopython/biopython/commit/f0f4536119947e7d4df838adf6283e545e0dee54) commit. So for now, you must either use biopython from git, patch your own version, or use the python2 version.
* This bug will cause the downloading of sequences to fail with an urllib error sometimes.

##Found a bug?
Or would like a feature added? Or maybe just wanto to drop some feedback?
Just open an issue on github!