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

import sys
import re
import os

from PyQt5 import QtGui, QtCore, QtWidgets

from back_end import Downloader
from sanity_checks import basic_checks


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.initUI()

    def initUI(self):
        ## Create widgets
        # Download button
        self.dlbtn = QtWidgets.QPushButton("Download!", self)
        self.dlbtn.setToolTip('Click to start downloading...')
        self.dlbtn.clicked.connect(self.statusChange)
        self.dlbtn.clicked.connect(self.runOnClick)
        self.dlbtn.resize(self.dlbtn.sizeHint())

        # Quit button
        self.qbtn = QtWidgets.QPushButton("Quit", self)
        self.qbtn.setToolTip('Click to exit the program...')
        self.qbtn.clicked.connect(QtCore.QCoreApplication.instance().quit)
        self.qbtn.resize(self.qbtn.sizeHint())

        # Cancel button
        self.canbtn = QtWidgets.QPushButton("Stop", self)
        self.canbtn.setToolTip('Stop the download.')
        self.canbtn.setIcon(QtGui.QIcon("assets/stop.png"))
        self.canbtn.resize(self.canbtn.sizeHint())
        self.canbtn.setEnabled(False)

        # Progress bar
        self.progbar = QtWidgets.QProgressBar(self)

        # Status bar
        self.statusBar().showMessage('Ready')

        # Title lable
        self.title = QtWidgets.QLabel(self)
        self.title.setText("NCBI Mass Sequence Downloader")
        self.title.setFont(QtGui.QFont("Sans", 16, QtGui.QFont.Bold, True))

        # Help with queries
        self.query_help_url = "http://www.ncbi.nlm.nih.gov/books/NBK3837/#_EntrezHelp_Entrez_Searching_Options_"
        self.manual_url = "http://ncbi-mass-sequence-downloader.readthedocs.org/en/latest"
        self.help_label = QtWidgets.QLabel(self)
        self.help_label.setText("Read the <a href=" + self.manual_url +
                                ">manual</a>. Help with the <a href=" +
                                self.query_help_url + ">query</a>.")
        self.help_label.setOpenExternalLinks(True)

        # Databases to search and respective label
        self.databases = QtWidgets.QComboBox(self)
        self.databases.addItem("nucleotide")
        self.databases.addItem("nuccore")
        self.databases.addItem("nucgss")
        self.databases.addItem("nucest")
        self.databases.addItem("protein")
        self.databases.addItem("genome")
        self.databases.addItem("popset")
        self.databases_label = QtWidgets.QLabel(self)
        self.databases_label.setText("Database to search:")

        # Search query line and respective label
        self.search_query = QtWidgets.QLineEdit(self)
        self.search_query.setFixedWidth(470)
        self.search_query_label = QtWidgets.QLabel(self)
        self.search_query_label.setText("Search Query:")

        # File management
        self.save_file_label = QtWidgets.QLabel(self)
        self.save_file_label.setText("File Location:")

        self.save_file_line = QtWidgets.QLineEdit(self)
        self.save_file_line.setFixedWidth(300)
        self.save_file_line.setEnabled(False)

        self.save_file_button = QtWidgets.QPushButton("Save as...", self)
        self.save_file_button.resize(self.dlbtn.sizeHint())
        self.save_file_button.setToolTip('Click to select file location...')
        self.save_file_button.clicked.connect(self.fileHandle)

        # Container widget
        self.main_widget = QtWidgets.QWidget(self)

        ## Set layout
        # Containers
        self.main_layout = QtWidgets.QVBoxLayout(self.main_widget)
        self.main_layout.sizeConstraint = QtWidgets.QLayout.SetDefaultConstraint

        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)

        # Box for progress bar
        self.progressBox = QtWidgets.QHBoxLayout()
        self.progressBox.addWidget(self.progbar)

        # Box for search query
        self.queryBox = QtWidgets.QHBoxLayout()
        self.queryBox.addWidget(self.search_query_label)
        self.queryBox.addWidget(self.search_query)
        self.queryBox.addStretch(1)

        # Box for query help and database
        self.query_database_box = QtWidgets.QHBoxLayout()
        self.query_database_box.addWidget(self.help_label)

        self.query_database_box.addStretch(1)

        self.query_database_box.addWidget(self.databases_label)
        self.query_database_box.addWidget(self.databases)

        # Box for file management
        self.file_box = QtWidgets.QHBoxLayout()
        self.file_box.addWidget(self.save_file_label)
        self.file_box.addWidget(self.save_file_line)
        self.file_box.addWidget(self.save_file_button)
        self.file_box.addWidget(self.canbtn)

        # Box for Title
        self.titlebox = QtWidgets.QHBoxLayout()
        self.titlebox.addStretch(1)
        self.titlebox.addWidget(self.title)
        self.titlebox.addStretch(1)

        # Box for bottom buttons
        self.bottomBox = QtWidgets.QHBoxLayout()
        self.bottomBox.addStretch(1)
        self.bottomBox.addWidget(self.dlbtn)
        self.bottomBox.addWidget(self.qbtn)

        # Vertical stack
        self.main_layout.addLayout(self.titlebox)
        self.main_layout.addStretch(1)
        self.main_layout.addLayout(self.query_database_box)
        self.main_layout.addLayout(self.queryBox)
        self.main_layout.addLayout(self.file_box)
        self.main_layout.addStretch(1)
        self.main_layout.addLayout(self.progressBox)
        self.main_layout.addLayout(self.bottomBox)

        # MainWindow proprieties
        self.setWindowTitle('NCBI mass downloader')
        self.setWindowIcon(QtGui.QIcon('assets/Icon.png'))

        # Draw it!
        self.show()


    def statusChange(self):
        if self.statusBar().currentMessage() == "Ready":
            self.statusBar().showMessage("Downloading...")
            self.canbtn.setEnabled(True)
            self.dlbtn.setEnabled(False)
            self.qbtn.setEnabled(False)
            self.save_file_button.setEnabled(False)
        elif self.statusBar().currentMessage() == "Downloading...":
            self.statusBar().showMessage("Ready")
            self.canbtn.setEnabled(False)
            self.dlbtn.setEnabled(True)
            self.qbtn.setEnabled(True)
            self.save_file_button.setEnabled(True)


    def fileHandle(self):
        self.savefile = QtWidgets.QFileDialog.getSaveFileName(self, "Save to file...", "", "Fasta Files (*.fasta);;All Files (*)")[0]
        self.save_file_line.setText(self.savefile)


    def runOnClick(self):
        # str() method used for python2 compatibility (python2 does not handle QString natively)
        self.database_to_search = str(self.databases.currentText())
        self.search_term = str(self.search_query.displayText())
        self.file_to_handle = str(self.save_file_line.displayText())


        if basic_checks_gui.sanity_checker(self.search_term, self.file_to_handle, 1) == 1:

            self.Get_data = DownloaderGui(self.database_to_search,
                                          self.search_term,
                                          self.file_to_handle,
                                          1)
            self.work_thread = QtCore.QThread()
            self.Get_data.max_seq.connect(self.progbar.setMaximum)
            self.Get_data.prog_data.connect(self.progbar.setValue)
            self.Get_data.no_match.connect(lambda msg: setattr(self, 'message', msg))
            self.Get_data.finished.connect(self.what_next)
            self.Get_data.moveToThread(self.work_thread)
            self.work_thread.started.connect(self.Get_data.run_everything)
            self.canbtn.clicked.connect(self.stop_threads)
            self.work_thread.start()


    def stop_threads(self):
        self.dl_message = "Download canceled."
        self.Get_data.terminated = True
        self.work_thread.terminate()
        self.work_thread.quit()
        self.statusBar().showMessage("Canceling the download. The program may become irresponsive for a while.")
        self.work_thread.wait()
        self.statusBar().showMessage("Downloading...")
        self.what_next(self.dl_message)


    def what_next(self, dl_message):
        if self.DlFinished(dl_message) == 2097152:
            self.close()
        else:
            self.cleanForms()
            self.statusChange()


    def cleanForms(self):
        """
        Clear forms for making another download
        """
        self.search_query.setText("")
        self.save_file_line.setText("")
        self.progbar.setValue(0)


    def DlFinished(self, message):
        """
        Create message box for finished download
        """
        self.question = QtWidgets.QMessageBox(self)
        self.question.setIcon(QtWidgets.QMessageBox.Question)
        self.question.setText(message)
        self.question.setInformativeText("Would you like to reset the forms and make another download or close the program?")
        self.question.setStandardButtons(QtWidgets.QMessageBox.Reset | QtWidgets.QMessageBox.Close)

        reply = self.question.exec_()

        return reply

    def sanityCheck(self):
        """
        Check if the variables to be sent to the back end make sense
        """
        if len(self.search_term) < 3:
            self.fail = QtWidgets.QMessageBox.warning(self, "Problem with search query", "Your search query is too short. It should have at least 3 characters.", QtWidgets.QMessageBox.Ok)
            return 0
        elif (os.path.exists(os.path.dirname(self.file_to_handle)) == False) or (os.access(os.path.dirname(self.file_to_handle), os.W_OK) == False):
            self.fail = QtWidgets.QMessageBox.warning(self, "Problem with path or permissions", "The path leading to your output file does not seem to exist or you don't have write permissions for it. Please correct it and try again.", QtWidgets.QMessageBox.Ok)
            return 0
        else:
            return 1


class DownloaderGui(Downloader, QtCore.QObject):
    # Just add PyQt 'magic' to Downloader() and create emmiters in constructor.
    prog_data = QtCore.pyqtSignal(int)
    max_seq = QtCore.pyqtSignal(int)
    no_match = QtCore.pyqtSignal(str)
    finished = QtCore.pyqtSignal(str)


    def __init__(self, database, term, outfile, gui):
        # Add threading!
        Downloader.__init__(self, database, term, outfile, gui)

class basic_checks_gui(basic_checks, QtCore.QObject):
    length_ok = QtCore.pyqtSignal(tuple)
    outfile_ok = QtCore.pyqtSignal(tuple)

    def __init__(self, query, outfile, gui):
        basic_checks.__inti__(self, query, outfile, gui)

def main():

    app = QtWidgets.QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
