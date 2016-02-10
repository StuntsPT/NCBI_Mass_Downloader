# Binary Builds

## Where to find them
Binary builds can always be found in the [releases page](https://github.com/StuntsPT/NCBI_Mass_Downloader/releases).
Note that not all versions get a binary build. Sometimes this happens because two versions are released too close to each other, other times, because the changes were only in the documentation.

## How to build them
Binaries for all platforms (GNU/Linux, Windows & OSX) are built using [pyinstaller](http://www.pyinstaller.org/). They are built using the following options:

(GNU/Linux)
```bash
pyinstaller -F -w NCBI_downloader.py
```

(Windows)
```powershell
C:\Python34\scripts\pyinstaller.exe -F -w -i assets/Icon.ico NCBI_downloader.py
```

(OSX)
```bash
pyinstaller -F -w -i assets/Icon.png NCBI_downloader.py
```

* The "-F" option will bundle all the used libs in a single file.
* The "-w" option will disable opening a console window to display STDOUT and STDERR (eg. the program will only spawn the GUI window).
* The "-i" option used in the Windows and OSX builds will bundle the picture in `assets` with the executable file.

The produced binaries are then bundled with the XSDs and DTDs directories, license files and README.md in a compressed file which is then uploaded to the respective release page in github.

For OSX, instead of bundling the executable file with the above mentioned files, the ".app" directory is bundled. This is done so that it's possible to run the program with a "double click" on OSX.

[Return to Introduction](index.md)
