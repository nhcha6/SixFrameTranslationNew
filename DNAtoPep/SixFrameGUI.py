from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, \
    QFileDialog, QGridLayout, QLabel, QComboBox, QCheckBox, QMessageBox, QDesktopWidget, \
    QProgressBar, QLineEdit, QInputDialog, QGroupBox, QFormLayout
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import *
from PyQt5.QtCore import pyqtSlot
import sys
from SixFrameTranslation import *
from time import time
import platform

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

        self.outputPath = None

    def closeEvent(self, event):
        print('closed')
        # windows close command
        if platform.system() == 'Windows':
            os.system('taskkill /f /fi "WINDOWTITLE eq Peptide Splicer" /t')
        # mac close command
        else:
            os.system("ps aux |grep MersGUI | grep -v 'pattern_of_process_you_dont_want_to_kill' | awk '{print $2}' |xargs kill")


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
        Called after generate output is clicked and the user confirms their input. Opens a window to select a file location
        to save the output to, and if valid opens a window to input the file name.
        """
        # opens a window to select file location.
        self.outputPath = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        # if no outout path is returned, simply return to the main GUI and the user can choose to recommence the file location
        # selection process if they desire.
        if self.outputPath == '':
            return
        # else if a valid path is selected, bring up a dialog to input the file name
        else:
            self.filePathDialog()

    def filePathDialog(self):
        """
        This function initialises and shows the filing naming popup.
        """
        self.outputNameBox = QGroupBox('Output Name')
        self.outputNameLayout = QFormLayout()
        self.outputNameLayout.addRow(QLabel("Add a name for the output file."))
        self.outputNameLayout.addRow(QLabel('Banned characters: \ / : * " < > |'))
        self.fileName = QLineEdit()
        self.fileName.textChanged[str].connect(self.nameChecker)
        self.button = QPushButton("Create Output")
        self.valid = QLabel("Valid")
        self.button.clicked.connect(self.returnPath)
        self.outputNameLayout.addRow(self.fileName, self.valid)
        self.outputNameLayout.addRow(self.button)
        self.outputNameBox.setLayout(self.outputNameLayout)
        self.outputNameBox.show()

    def nameChecker(self, input):
        """
        This function is called every time the file name lineEdit is updated. It takes the param input, which is the
        text in the lineEdit, and checks if it is a valid file name.
        :param input:
        """
        # assign bannedCharacters to variables.
        bannedCharacters = set('\/:*"<>|')
        # if the input has no intersection with the banned characters it is valid. If so, update the label validity label
        # and set ensure the generate output button is concurrently enabled/disabled.
        if len(set(input).intersection(bannedCharacters)) == 0:
            self.valid.setText("Valid")
            self.button.setEnabled(True)
        else:
            self.valid.setText("Invalid")
            self.button.setEnabled(False)

    def returnPath(self):
        start = time()
        # create the output file name by combining path with the name.
        outputFile = self.outputPath + '/' + self.fileName.text()
        print(outputFile)
        self.outputGen = OutputGenerator(self.createOutput, outputFile, self.minPeptideLen, self.inputFile,
                                         self.removeSubFlag, self.writeSubFlag, self.originFlag)
        self.outputGen.signals.finished.connect(self.outputFinished)
        self.threadpool.start(self.outputGen)
        self.outputLabel = QLabel("Generating Output. Please Wait!")
        self.grid.addWidget(self.outputLabel, 7, 1)
        # generateOutputNew(outputPath, self.minPeptideLen, self.inputFile)
        end = time()
        print(end - start)
        # close the output name box.
        self.outputNameBox.close()

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
                self.getOutputPath()

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