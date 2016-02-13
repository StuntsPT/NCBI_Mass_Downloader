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

import os

class basic_checks(object):
    def __inti__(self, query, outfile, gui):
        self.query = query
        self.outfile = outfile
        self.gui = gui
        super(basic_checks, self).__init__()


    def sanity_checker(self):
        """
        Make some basic variable checks.
        """
        if len(self.query) < 3:
            self.msg = "Your search query is too short. It should have at least 3 "\
                  "characters."
            if self.gui == 1:
                self.length_ok.emit("Problem with search query", msg)
                return 0
            else:
                quit(self.msg)

        elif (os.path.exists(os.path.dirname(self.outfile)) == False) or (os.access(os.path.dirname(self.outfile), os.W_OK) == False):
            self.msg = "The path leading to your output file does not seem to exist or "\
                  "you don't have write permissions for it. Please correct it and "\
                  "try again."
            if self.gui == 1:
                self.fail = QtWidgets.QMessageBox.warning(self, "Problem with path or permissions", self.msg, QtWidgets.QMessageBox.Ok)
                return 0
            else:
                quit(self.msg)

        return 1
