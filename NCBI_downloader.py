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
import re

from PyQt4 import QtGui, QtCore
from Bio import Entrez
from os import remove, stat
from shutil import move

class MainWindow(QtGui.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.initUI()

    def initUI(self):
        ##Create widgets
        #Download button
        dlbtn = QtGui.QPushButton("Download!", self)
        dlbtn.setToolTip('Click to start downloading...')
        dlbtn.clicked.connect(self.statusChange)
        dlbtn.clicked.connect(self.runOnClick)
        dlbtn.resize(dlbtn.sizeHint())

        #Quit button
        qbtn = QtGui.QPushButton("Quit", self)
        qbtn.setToolTip('Click exit the program...')
        qbtn.clicked.connect(QtCore.QCoreApplication.instance().quit)
        qbtn.resize(qbtn.sizeHint())

        #Progress bar
        self.progbar = QtGui.QProgressBar(self)
        #self.pbar.setValue(self."X") Reminder!

        #Status bar
        self.statusBar().showMessage('Ready')

        #Title lable
        self.title = QtGui.QLabel(self)
        self.title.setText("NCBI mass sequence downloader")
        self.title.setFont(QtGui.QFont("Sans", 16, QtGui.QFont.Bold, True))

        #Email line edit and respective label
        self.email_line = QtGui.QLineEdit(self)
        self.email_line.setFixedWidth(220)
        self.email_label = QtGui.QLabel(self)
        self.email_label.setText("Email address:")

        #Databases to search and respective label
        self.databases = QtGui.QComboBox(self)
        self.databases.addItem("nucleotide")
        self.databases.addItem("nuccore")
        self.databases.addItem("nucgss")
        self.databases.addItem("nucest")
        self.databases.addItem("protein")
        self.databases.addItem("genome")
        self.databases.addItem("popset")
        self.databases_label = QtGui.QLabel(self)
        self.databases_label.setText("Database to search:")

        #Search query line and respective label
        self.search_query = QtGui.QLineEdit(self)
        self.search_query.setFixedWidth(300)
        self.search_query_label = QtGui.QLabel(self)
        self.search_query_label.setText("Search Query:")

        #File management
        self.save_file_label = QtGui.QLabel(self)
        self.save_file_label.setText("File Location:")

        self.save_file_line = QtGui.QLineEdit(self)
        self.save_file_line.setFixedWidth(300)

        self.save_file_button = QtGui.QPushButton("Find...", self)
        self.save_file_button.resize(dlbtn.sizeHint())
        self.save_file_button.setToolTip('Click to select file location...')
        self.save_file_button.clicked.connect(self.fileHandle)

        #self.savefile = QtGui.QFileDialog.getSaveFileName(self, "Save to file...", "", ".fasta")

        #Container widget
        self.main_widget = QtGui.QWidget(self)

        ##Set layout
        #Containers
        self.main_layout = QtGui.QVBoxLayout(self.main_widget)
        self.main_layout.sizeConstraint = QtGui.QLayout.SetDefaultConstraint

        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)

        #Box for progress bar
        self.progressBox = QtGui.QHBoxLayout()
        self.progressBox.addWidget(self.progbar)

        #Box for search query
        self.queryBox = QtGui.QHBoxLayout()
        self.queryBox.addWidget(self.search_query_label)
        self.queryBox.addWidget(self.search_query)
        self.queryBox.addStretch(1)

        #Box for email and database
        self.email_database_box = QtGui.QHBoxLayout()
        self.email_database_box.addWidget(self.email_label)
        self.email_database_box.addWidget(self.email_line)

        self.email_database_box.addStretch(1)

        self.email_database_box.addWidget(self.databases_label)
        self.email_database_box.addWidget(self.databases)

        #Box for file management
        self.file_box = QtGui.QHBoxLayout()
        self.file_box.addWidget(self.save_file_label)
        self.file_box.addWidget(self.save_file_line)
        self.file_box.addWidget(self.save_file_button)

        #Box for Title
        self.titlebox = QtGui.QHBoxLayout()
        self.titlebox.addStretch(1)
        self.titlebox.addWidget(self.title)
        self.titlebox.addStretch(1)

        #Box for bottom buttons
        self.bottomBox = QtGui.QHBoxLayout()
        self.bottomBox.addStretch(1)
        self.bottomBox.addWidget(dlbtn)
        self.bottomBox.addWidget(qbtn)

        #Vertical stack
        self.main_layout.addLayout(self.titlebox)
        self.main_layout.addStretch(1)
        self.main_layout.addLayout(self.email_database_box)
        self.main_layout.addLayout(self.queryBox)
        self.main_layout.addLayout(self.file_box)
        self.main_layout.addStretch(1)
        self.main_layout.addLayout(self.progressBox)
        self.main_layout.addLayout(self.bottomBox)


        #MainWindow proprieties
        self.setWindowTitle('NCBI mass downloader')
        self.setWindowIcon(QtGui.QIcon('assets/Icon.png'))

        #Draw it!
        self.show()

    def closeEvent(self, event):

        reply = QtGui.QMessageBox.question(self, 'Really quit?',
            "Are you sure you want to quit?", QtGui.QMessageBox.Yes |
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def statusChange(self):
        if self.statusBar().currentMessage() == "Ready":
            self.statusBar().showMessage("Downloading...")
        elif self.statusBar().currentMessage() == "Downloading...":
            self.statusBar().showMessage("Ready")

    def fileHandle(self):
        self.savefile = QtGui.QFileDialog.getSaveFileName(self, "Save to file...", "", "Fasta Files (*.fasta);;All Files (*)")
        self.save_file_line.setText(self.savefile)

    def runOnClick(self):
        self.email_address = self.email_line.displayText()
        self.database_to_search = self.databases.currentText()
        self.search_term = self.search_query.displayText()
        self.file_to_handle = self.save_file_line.displayText()

        #TODO: Implement argument detection. (trust no-one)

        Get_data = DownloaderGui(self.email_address, self.database_to_search, self.search_term, self.file_to_handle, 1)
        Get_data.max_seq.connect(self.progbar.setMaximum)
        Get_data.prog_data.connect(self.progbar.setValue)
        Get_data.runEverything()
        #return self.email_address, self.database_to_search, self.search_term, self.file_to_handle, 1

class Downloader():
        def __init__(self, email, database, term, outfile, gui):
            self.email = email
            self.database = database
            self. term = term
            self.outfile = outfile
            self.gui = gui
            #super(Downloader, self).__init__()

        def NCBI_search(self):
            #Submit search to NCBI and return the records
            handle = Entrez.esearch(db=self.database,term=self.term,usehistory="y",retmax=10000000)
            record = Entrez.read(handle)
            handle.close()

            return record

        def NCBI_post(self, IDs):
            #Submit id_list to NCBI via epost and return the records
            IDs_string = ",".join(IDs)
            handle = Entrez.epost(self.database,id=IDs_string,retmax=10000000)
            record = Entrez.read(handle)
            handle.close()

            count = len(IDs)
            webenv = record["WebEnv"]
            query_key = record["QueryKey"]

            return count, IDs, webenv, query_key

        def Record_processor(self,record):
            #Processes the record into sparate usefull information
            count = int(record["Count"])
            IDs = record["IdList"]
            webenv = record["WebEnv"]
            query_key = record["QueryKey"]

            assert count == len(IDs)

            return count, IDs, webenv, query_key

        def NCBI_history_fetch(self, count, IDs, webenv, query_key, Bsize, Run):
            #Fetches results from NCBI using history
            try:
                a = open(self.outfile,'r')
                a.close()
            except:
                a = open(self.outfile,'w')
                a.close()
            if Run == 1 and stat(self.outfile).st_size != 0:
                self.ReDownloader(IDs)
            else:

                outfile = open(self.outfile,'a')
                if self.gui == 1:
                        self.max_seq.emit(count)
                if Bsize > count:
                    Bsize = count
                for start in range(0,count,Bsize):
                    if start + Bsize < count:
                        end = start + Bsize
                    else:
                        end = count
                    print("Downloading record %i to %i of %i" %(start+1, end, count))

                    if self.gui == 1:
                        self.prog_data.emit(start)

                    #Make sure that even on server errors the program carries on.
                    #If the servers are dead, well, you were not going anywhere anyway...
                    while True:
                        try:
                            fetch_handle = Entrez.efetch(db=self.database, rettype="fasta", retstart=start, retmax=Bsize, webenv=webenv, query_key=query_key)
                            break
                        except:
                            pass
                    data = fetch_handle.read()
                    fetch_handle.close()
                    outfile.write(data)

            outfile.close()
            self.ReDownloader(IDs)

        def ReDownloader(self, IDs):
            #Check for missing sequences:
            print("Checking for sequences that did not download... Please wait.")
            ver_IDs = self.Error_finder()
            missing_IDs = set()
            for i in IDs:
                if i not in ver_IDs:
                    missing_IDs.add(i)
            IDs = missing_IDs #Improve performance on subsequent runs
            if len(missing_IDs) == 0:
                print("All sequences were downloaded correctly. Good!")
                if self.gui == 0:
                    quit("Program finished without error.")
            else:
                print("%s sequences did not download correctly (or at all). Retrying..." %(len(missing_IDs)))
                count, IDs, webenv, query_key = self.NCBI_post(IDs)
                self.NCBI_history_fetch(count, IDs, webenv, query_key, 1000, 2)

        def Error_finder(self):
            #Looks for errors in the output fasta and retruns a list of necessary retries
            temp_file = self.outfile + ".tmp"
            move(self.outfile, temp_file)
            original_file = open(temp_file,'r')
            new_file = open(self.outfile,'w')
            verified_IDs = set()

            for lines in original_file:
                if lines.startswith(">"):
                    ID = re.search("gi\|.*?\|",lines).group(0)[3:-1]
                    verified_IDs.add(ID)
                    new_file.write("\n" + lines) #TODO: remove first empty line from file
                elif lines.strip().startswith("<") or lines.startswith("\n"):
                    pass
                else:
                    new_file.write(lines)

            original_file.close()
            new_file.close()
            remove(temp_file)

            return verified_IDs

        def runEverything(self):
            #Run the functions in order
            batch_size = 1000
            Entrez.email = self.email

            rec = self.NCBI_search()
            count, IDs, webenv, query_key = self.Record_processor(rec)
            self.NCBI_history_fetch(count, IDs, webenv, query_key, batch_size, 1)

class DownloaderGui(Downloader, QtCore.QThread):
    prog_data = QtCore.pyqtSignal(int)
    max_seq = QtCore.pyqtSignal(int)





def main():

    if len(sys.argv) < 2:
        app = QtGui.QApplication(sys.argv)
        ex = MainWindow()
        sys.exit(app.exec_())
        print("OK")
    else:
        dl = Downloader(sys.argv[1],sys.argv[3],sys.argv[2],sys.argv[4], 0)
        dl.runEverything()


if __name__ == '__main__':
    main()