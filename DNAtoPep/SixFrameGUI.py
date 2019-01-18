from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, \
    QFileDialog, QGridLayout, QLabel, QComboBox, QCheckBox, QMessageBox, QDesktopWidget, \
    QProgressBar, QLineEdit, QInputDialog
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import *
from PyQt5.QtCore import pyqtSlot
import sys
from SixFrameTranslation import *
from time import time

class WorkerSignals(QObject):
    """
    Signals class that is used for the GUI when emitting custom signals
    """

    finished = pyqtSignal()

class OutputGenerator(QRunnable):
    """

    """

    def __init__(self, fn, *args, **kwargs):
        super(OutputGenerator, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        self.fn(*self.args)
        self.signals.finished.emit()


class Example(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()
        self.inputFile = ""
        self.minPeptideLen = 0
        self.originFlag = ""
        self.writeSubFlag = ""
        self.removeSubFlag = ""
        self.threadpool = QThreadPool()

    def initUI(self):

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.initialiseWidgets()

        self.move(300, 150)
        self.setWindowTitle('DNA to Protein')
        self.show()

    def initialiseWidgets(self):
        self.importDNA = QPushButton('Import DNA fasta')
        self.grid.addWidget(self.importDNA, 1, 1)
        self.importDNA.clicked.connect(self.uploadInput)

        self.minLenCombo = QComboBox()
        self.minLenLabel = QLabel('Minimum Protein Length: ')
        self.grid.addWidget(self.minLenCombo, 2,2)
        self.grid.addWidget(self.minLenLabel, 2,1)
        for i in range(2, 15):
             self.minLenCombo.addItem(str(i))

        self.removeSubseq = QCheckBox('Remove Subsequences')
        self.removeSubseq.stateChanged.connect(self.disableCheckboxes)
        self.ignoreOrigin = QCheckBox('Ignore Origin')
        self.ignoreOrigin.setEnabled(False)
        self.writeSubseq = QCheckBox('Print Deleted Subseq')
        self.writeSubseq.setEnabled(False)
        self.grid.addWidget(self.removeSubseq, 3, 1)
        self.grid.addWidget(self.ignoreOrigin, 4, 1)
        self.grid.addWidget(self.writeSubseq, 3, 2)

        self.generateOutput = QPushButton('Generate Output')
        self.grid.addWidget(self.generateOutput, 5,1)
        self.generateOutput.clicked.connect(self.outputCheck)

    def uploadInput(self):
        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home/')
        if fname[0][-5:] == 'fasta':
            print('in')
            self.inputFile = fname[0]
            QMessageBox.about(self, 'Message', 'Fasta input file successfully uploaded!')

    def getOutputPath(self):

        """
        Called after generate output is clicked. Opens a window to select a file location to save the output to.
        Returns False if no path is selected, otherwise returns the selected path.
        """

        outputFile = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        if outputFile == '':
            return False
        else:
            text, ok = QInputDialog.getText(self, 'Input Dialog',
                                            'Enter your file name:')

            if ok:
                outputPath = outputFile + '/' + text
            else:
                return False
        print(outputPath)
        return outputPath

    def outputCheck(self):
        if self.inputFile == "":
            QMessageBox.about(self, 'Message', 'Please Upload a Fasta File before generating output!')
        else:
            minString= self.minLenCombo.currentText()
            self.minPeptideLen = int(minString)
            self.removeSubFlag = self.removeSubseq.isChecked()
            self.writeSubFlag = self.writeSubseq.isChecked()
            self.originFlag = self.ignoreOrigin.isChecked()
            reply = QMessageBox.question(self, 'Message', 'Do you wish to confirm the following input?\n' +
                                          'Minimum Protein Length: ' + minString + '\n' +
                                          'Remove Subsequences: ' + str(self.removeSubFlag) + '\n' +
                                          'Write Subsequences to File: ' + str(self.writeSubFlag) + '\n' +
                                          'Ignore Origin Sequences: ' + str(self.originFlag) + '\n' +
                                          'Input File ' + self.inputFile)
            if reply == QMessageBox.Yes:
                start = time()
                outputPath = self.getOutputPath()
                if outputPath is not False:

                    self.outputGen = OutputGenerator(self.createOutput, outputPath, self.minPeptideLen, self.inputFile, self.removeSubFlag, self.writeSubFlag, self.originFlag)
                    self.outputGen.signals.finished.connect(self.outputFinished)
                    self.threadpool.start(self.outputGen)
                    self.outputLabel = QLabel("Generating Output. Please Wait!")
                    self.grid.addWidget(self.outputLabel,7,1)
                    #generateOutputNew(outputPath, self.minPeptideLen, self.inputFile)
                    end = time()
                    print(end-start)

    def createOutput(self, outputPath, minPeptideLen, inputFile, removeSubFlag, writeSubFlag, originFlag):
        generateOutputNew(outputPath, self.minPeptideLen, self.inputFile, removeSubFlag, writeSubFlag, originFlag)

    def outputFinished(self):
        QMessageBox.about(self, "Message", "All done!")
        self.grid.removeWidget(self.outputLabel)
        self.outputLabel.deleteLater()
        self.outputLabel = None

    def disableCheckboxes(self):
        if self.removeSubseq.isChecked():
            self.writeSubseq.setEnabled(True)
            self.ignoreOrigin.setEnabled(True)
        else:
            self.writeSubseq.setChecked(False)
            self.ignoreOrigin.setChecked(False)
            self.writeSubseq.setEnabled(False)
            self.ignoreOrigin.setEnabled(False)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())