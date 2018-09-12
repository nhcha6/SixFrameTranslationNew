from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, \
    QFileDialog, QGridLayout, QLabel, QComboBox, QCheckBox, QMessageBox, QDesktopWidget, \
    QProgressBar, QLineEdit
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import *
from PyQt5.QtCore import pyqtSlot
import sys
#from 6FrameTranslation.py import *


class Example(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        # names = ['Cls', 'Bck', '', 'Close',
        #          '7', '8', '9', '/',
        #          '4', '5', '6', '*',
        #          '1', '2', '3', '-',
        #          '0', '.', '=', '+']
        #
        # positions = [(i, j) for i in range(5) for j in range(4)]
        #
        # for position, name in zip(positions, names):
        #
        #     if name == '':
        #         continue
        #     button = QPushButton(name)
        #     grid.addWidget(button, *position)

        self.initialiseWidgets()

        self.move(300, 150)
        self.setWindowTitle('DNA to Protein')
        self.show()

    def initialiseWidgets(self):
        self.importDNA = QPushButton('Import DNA fasta')
        self.grid.addWidget(self.importDNA, 1, 1)
        self.importDNA.clicked.connect(self.uploadInput)

    def uploadInput(self):
        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home/')
        print(fname)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())