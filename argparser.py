#!/usr/bin/python3
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

import argparse
# Argument list
def arg_list():
    parser = argparse.ArgumentParser(description="A program to download large "
                                                 "amounts of sequences from "
                                                 "NCBI databases. Find the "
                                                 "complete manual here: "
                                                 "http://ncbi-mass-sequence-downloader.readthedocs.org/en/latest/",
                                     prog="NCBI_downloader.py",
                                     formatter_class=argparse.RawTextHelpFormatter)

    io_opts = parser.add_argument_group("Input/Output options")
    query_opts = parser.add_argument_group("Query options")


    io_opts.add_argument("-o", "--out",
                         dest="outfile",
                         required=True,
                         type=str,
                         default=None,
                         metavar="filepath",
                         help="Path to where you wish to save your sequences.")

    query_opts.add_argument("-d", "--db",
                            dest="database",
                            required=True,
                            type=str,
                            default=None,
                            choices=["nucleotide", "nuccore", "nucgss",
                                     "nucest", "protein", "genome", "popset"],
                            metavar="database",
                            help='Database to query. Valid values are '
                                 '"nucleotide", "nuccore", "nucgss", "nucest", '
                                 '"protein", "genome" and "popset".')

    query_opts.add_argument("-q", "--query",
                            dest="query",
                            required=True,
                            type=str,
                            default=None,
                            metavar='"query term"',
                            help="Your query term. For more information on how "
                                 "to query the Entrez database, please look "
                                 "here: http://www.ncbi.nlm.nih.gov/books/NBK3837/#_EntrezHelp_Entrez_Searching_Options_ ."
                            )

    arg = parser.parse_args()
    if len(arg.query) < 3:
        quit_text = "usage: NCBI_downloader.py [-h] -o filepath -d database "\
                    "-q query term\nNCBI_downloader.py: error: argument "\
                    "-q/--query: invalid query: '" + arg.query + "' (query "\
                    "must be at least 3 characters long)"
        quit(quit_text)

    return (arg.database, arg.query, arg.outfile)
