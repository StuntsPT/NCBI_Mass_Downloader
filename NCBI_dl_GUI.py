#!/usr/bin/python2

import sys
from PyQt4 import QtGui, QtCore

class MainWindow(QtGui.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.initUI()

    def initUI(self):

        dlbtn = QtGui.QPushButton('Download!', self)
        dlbtn.setToolTip('Click to start downloading...')
        dlbtn.clicked.connect(self.statusChange)
        dlbtn.resize(dlbtn.sizeHint())
        dlbtn.move(50, 50)

        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.setToolTip('Click exit the program...')
        qbtn.clicked.connect(QtCore.QCoreApplication.instance().quit)
        qbtn.resize(qbtn.sizeHint())
        qbtn.move(50, 100)

        self.statusBar().showMessage('Ready')

        self.resize(250, 150)
        self.setWindowTitle('NCBI mass downloader')
        self.setWindowIcon(QtGui.QIcon('assets/Icon.png'))

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
def main():

    app = QtGui.QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())



if __name__ == '__main__':
    main()