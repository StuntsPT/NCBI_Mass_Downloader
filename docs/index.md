# NCBI Mass Sequence Downloader

## Description
This program will download sequences *en masse* from several NCBI databases (at the user's chioce).
After the downloading is finished, the program will check the resulting file for any missing sequences and continuously retry the download until all sequences are present in the local file.
It will run on both python2 and python3.
The program uses the Entrez module from [biopython](https://github.com/biopython/biopython). It is slightly modified from the original it was forked from. A big Kudos to the authors.


[NCBI Data Usage Policies and Disclaimers](https://www.ncbi.nlm.nih.gov/home/about/policies.shtml) may apply to downloaded data.

## Dependencies
* Pyhton (2 or 3) (and the standard lib);
* [PyQt5](http://www.riverbankcomputing.com/software/pyqt/intro) for the GUI version.

## Where to get it
* Source code - [NCBI Mass Sequence Downloader on github](https://github.com/StuntsPT/NCBI_Mass_Downloader)
* [Binaries](https://github.com/StuntsPT/NCBI_Mass_Downloader/releases) for multiple platforms.

## Contents
* [Installation & dependencies](install.md)
* [Binary builds](binaries.md)
* [Usage](usage.md)
* [Testing](testing.md)
* [Future plans](future.md)
* [Limitations](limits.md)
* [FAQ](faq.md)

## Bug reporting
Found a br or would like a feature added? Or maybe drop some feedback?
Just [open a new issue](https://github.com/StuntsPT/NCBI_Mass_Downloader/issues/new).

## License
GPLv3

Biopython files (Entrez.py, Parser.py and py3k.py) are under the [Biopython License Agreement](http://www.biopython.org/DIST/LICENSE)
