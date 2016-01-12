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

#Usage: NCBI_downloader.py "user@email-address.com" "Query term" "database" outfile.fasta
#Or run without arguments for the GUI version.

import sys
import signal

def kill_switch(*args):
    """
    Make sure the program will allways exit on Ctrl+C.
    """
    print("Ctrl+C detected, ending all processes and quitting.")
    sys.exit(0)

def main():

    if len(sys.argv) < 2:
        import front_end
        front_end.main()

    else:
        from back_end import Downloader
        dler = Downloader(sys.argv[1], sys.argv[3], sys.argv[2], sys.argv[4], 0)
        dler.run_everything()

if __name__ == '__main__':
    signal.signal(signal.SIGINT, kill_switch)
    main()
