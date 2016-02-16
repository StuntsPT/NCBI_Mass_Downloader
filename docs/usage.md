# Usage
To use the GUI version:

    python NCBI_downloader.py

To use the command line version:

    python NCBI_downloader.py -q "Query term" -d "database" -o outfile.fasta

## Notes
* The program will not overwrite anything in the output file, but rather append sequences to it;s
    * If you do not want this behavior, please use an empty file as output.
* On the top left of the GUI there is a link to [documentation on how to perform a query to the Entrez services](http://www.ncbi.nlm.nih.gov/books/NBK3837/#_EntrezHelp_Entrez_Searching_Options_).
* You may use *NCBI Mass Seqence Downloader* to access any of these databases: "nucleotide", "nuccore", "nucgss", "nucest", "protein", "genome" and "popset".
* *Query term* can take any argument just like in the *Entrez* website search engine (eg. "Lacerta monticola[organism]");

[Return to Introduction](index.md)
