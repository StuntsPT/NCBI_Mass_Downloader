#!/usr/bin/python2

import sys
from PyQt4 import QtGui, QtCore

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
        self.savefile = QtGui.QFileDialog.getSaveFileName(self, "Save to file...", "", "*.fasta")
        return self.savefile


def main():

    app = QtGui.QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()