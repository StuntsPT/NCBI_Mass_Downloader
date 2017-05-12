# NCBI Mass Sequence Downloader


## Description:
This program will download sequences *en masse* from several NCBI databases (at the user's choice).
After the downloading is finished, the program will check the resulting file for any missing sequences and continuously retry the download until all sequences are present in the local file.

[NCBI Data Usage Policies and Disclaimers](https://www.ncbi.nlm.nih.gov/home/about/policies.shtml) may apply to downloaded data.

## Documentation:
[![Documentation Status](https://readthedocs.org/projects/ncbi-mass-sequence-downloader/badge/?version=latest)](http://ncbi-mass-sequence-downloader.readthedocs.io/en/latest/?badge=latest)

The complete manual can be found on [readthedocs.org](http://ncbi-mass-sequence-downloader.readthedocs.io/en/latest/).

## Citation:
If you use NCBI Mass Sequence Downloader in your research, please cite the following paper:
Pina-Martins, F., & Paulo, O. S. (2016). NCBI Mass Sequence Downloader–Large dataset downloading made easy. SoftwareX, 5, 80–83. https://doi.org/10.1016/j.softx.2016.04.007


## Dependencies:
* Pyhton (2 or 3) (and the standard lib);
* [PyQt5](http://www.riverbankcomputing.com/software/pyqt/intro) for the GUI version.


## License:
GPLv3

Biopython files (Entrez.py, Parser.py and py3k.py) are under the [Biopython License Agreement](http://www.biopython.org/DIST/LICENSE)

## Found a bug?
Or would like a feature added? Or maybe drop some feedback?
Just [open a new issue](https://github.com/StuntsPT/NCBI_Mass_Downloader/issues/new).
