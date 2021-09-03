import os
import sys
import json
import re
import time
import shutil
import threading
import requests
import PyQt5
from PyQt5 import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication
import bs4
from Bio import Entrez
from Bio import SeqIO
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
import tpl_search_temp
import run_pfamscan
import run_ncbi_blast
import run_fasta_blast
import run_psi_blast
import run_psi_search


class TempSearch(QMainWindow, tpl_search_temp.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    params = {}
    search_result = {}
    aaLettersL = []
    accessionL = []
    expectL = [1e-300, 1e-100, 1e-50, 1e-10, 1e-5, 0.0001, 0.001, 0.1, 1, 2, 5, 10, 20, 50]
    blastTypeL = ['NCBI BLASTp', 'PSI-Search'] # 'FASTA BLAST', 'FASTM BLAST', 'PSI-BLAST',
    searchTypeL = ['Whole sequence', 'Residue range', 'Functional domain']
    blastParamL = ['NCBI BLASTp Parameters', 'PSI-Search Parameters'] # 'FASTA BLAST Parameters', 'FASTM BLAST Parameters', 'PSI-BLAST Parameters',
    blastSelectedL = []
    resRangeL = []
    domainDic = {}

    # _______ ncbi blastp _______

    params_ncbi = {}
    expectL_NCBI = [1e-200, 1e-100, 1e-50, 1e-10, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10, 100, 1000]
    matrixL_NCBI = ['BLOSUM62', 'PAM30', 'PAM70', 'PAM250', 'BLOSUM80', 'BLOSUM45', 'BLOSUM50', 'BLOSUM90']
    gapOpenBL62L = [6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 19]
    gapOpenPAM30L = [5, 6, 7, 8, 9, 10]
    gapOpenPAM70L = [6, 7, 8, 9, 10, 11]
    gapOpenPAM250L = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    gapOpenBL80L = [6, 7, 8, 9, 10, 11, 13, 25]
    gapOpenBL45L = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    gapOpenBL50L = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    gapOpenBL90L = [5, 6, 7, 8, 9, 10]
    gapExtL_NCBI = [0, 1, 2, 3, 4]
    compositionalL = ['F', 'D', 1, 2, 3]
    yesNoL = ['Yes', 'No']
    filterAbbrL = ['T', 'F']
    scoreL_NCBI = [5, 10, 20, 50, 100, 150, 200, 250, 500, 750, 1000]
    DropoffL = [0, 2, 4, 6, 8, 10]
    gapAlignL = ['true', 'false']

    # _______ fasta blast _______

    params_fasta = {}
    programL = ['FASTA', 'FASTX', 'FASTY', 'GGSEARCH', 'GLSEARCH']
    expectUpL = ['1e-600', 1e-300, 1e-100, 1e-50, 1e-10, 1e-5, 0.001, 0.1, 1.0, 2, 5, 10, 20, 50]
    expectL_FASTA = [0, 1e-300, 1e-100, 1e-50, 1e-10, 1e-5, 0.001, 0.1, 1.0, 2, 5, 10, 20, 50]
    matrixL_FASTA = ['BLOSUM50', 'BLASTP62', 'BLOSUM80', 'PAM250', 'PAM120', 'MDM40', 'MDM20', 'MDM10', 'VTML160',
                     'VTML120', 'VTML80', 'VTML40', 'VTML20', 'VTML10']
    matrixAbbrL_FASTA = ['BL50', 'BP62', 'BL80', 'P250', 'P120', 'M40', 'M20', 'M10', 'VT160', 'VT120', 'VT80',
                         'VT40', 'VT20', 'VT10']
    gapOpenL_FASTA = sorted([-30, -32, -35, -40, -50, -64] + [i for i in range(-25, 1)], reverse=True)
    gapExtL_FASTA = sorted([-16] + [i for i in range(-8, 1)], reverse=True)
    ktupL = [1, 2]
    #ktupAbbrL = [6, 5, 4, 3, 2, 1, -1]
    filterL_FASTA = ['None', 'seg', 'xnu', 'seg+xnu']
    scoreL_FASTA = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 500, 750, 1000]
    statisticL = ['Regress', 'MLE', 'Altshul-Gish', 'Regress/shuf.', 'MLE/shuf.']
    statisticAbbrL = [1, 2, 3, 11, 12]
    trFaL = ['true', 'false']

    # _______ fastm _______
    programL_fastm = ['FASTM', 'FASTF', 'FASTS']
    matrixL_FASTM = ['BLOSUM50', 'BLOSUM62', 'BLASTP62', 'BLOSUM80', 'PAM120', 'PAM250', 'MDM40', 'MDM20', 'MDM10']
    gapOpenL_FASTM = sorted([i for i in range(-23, 1)], reverse=True)
    gapExtL_FASTM = sorted([i for i in range(-8, 1)], reverse=True)

    # _______ psi blast _______

    params_psi_blast = {}
    pssmL = [1.0e-6, 1.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2, 2.0e-2, 0.1, 0.3, 0.5, 1.0, 3.0,
             10.0]
    expectL_PSIBL = [1.0e-200, 1.0e-100, 1.0e-50, 1.0e-10, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 0.1, 1.0, 10, 100, 1000]
    matrixL_PSIBL = ['BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'PAM30', 'PAM70']
    gapOpenBL45L_PSIBL = [10, 11, 12, 13, 14, 15, 16]
    gapOpenBL62L_PSIBL = [8, 9, 10, 11, 12, 13, 15, 16]
    gapOpenBL80L_PSIBL = [8, 9, 10, 11, 13]
    gapOpenPAM30L_PSIBL = [8, 9, 10]
    gapOpenPAM70L_PSIBL = [8, 9, 10, 11]
    gapExtL_PSIBL = [0, 1, 2, 3]
    filterAbbrL_PSIBL = ['T', 'F']
    scoreL_PSIBL = [5, 10, 20, 50, 100, 150, 200, 250, 500, 750, 1000, 5000]
    DropoffL_PSIBL = [0, 2, 4, 6, 8, 10, 15, 20, 25, 30]
    gapAlignL_PSIBL = ['true', 'false']
    cpFile = ''

    # _______ psi search _______

    params_psi_search = {}
    expectL_PSISE = [1.0e-200, 1.0e-100, 1.0e-50, 1.0e-10, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 0.1, 1.0, 10]
    matrixL_PSISE = ['BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'PAM30', 'PAM70']
    gapOpenBL45L_PSISE = [10, 11, 12, 13, 14, 15, 16]
    gapOpenBL62L_PSISE = [8, 9, 10, 11, 12, 13, 15, 16]
    gapOpenBL80L_PSISE = [8, 9, 10, 11, 13]
    gapOpenPAM30L_PSISE = [8, 9, 10]
    gapOpenPAM70L_PSISE = [8, 9, 10, 11]
    gapExtL_PSISE = [0, 1, 2, 3]
    trueFalse = ['true', 'false']
    filterL_PSISE = ['None', 'seg', 'xnu', 'seg+xnu']
    scoreL_PSISE = [5, 10, 20, 50, 100, 150, 200, 250, 500, 750, 1000, 5000]

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        # browse seq file
        self.browsSeq.released.connect(self.browse_seq)
        # blast type
        self.blastType.addItems(self.blastTypeL)
        # search type
        self.searchType.addItems(self.searchTypeL)
        self.searchType.currentTextChanged.connect(self.search_type)
        # pfamscan expect
        self.expect.addItems(str(i) for i in self.expectL)
        self.expect.setCurrentIndex(11)
        # blast parameter type
        self.blastParam.addItems(self.blastParamL)
        self.blastParam.currentTextChanged.connect(self.show_blast_param)
        # blast parameter list
        self.group_NCBI.show()
        self.group_FASTA.hide()
        self.horizontalLayoutWidget_40.setGeometry(QtCore.QRect(0, 0, 0, 0))
        self.group_FASTM.hide()
        self.horizontalLayoutWidget_92.setGeometry(QtCore.QRect(0, 0, 0, 0))
        self.group_PSI_BLAST.hide()
        self.horizontalLayoutWidget_69.setGeometry(QtCore.QRect(0, 0, 0, 0))
        self.group_PSI_search.hide()
        self.horizontalLayoutWidget_88.setGeometry(QtCore.QRect(0, 0, 0, 0))
        # ___ NCBI parameter value ___
        self.matrixTypesNCBI.addItems(self.matrixL_NCBI)
        self.matrixTypesNCBI.currentTextChanged.connect(self.get_gap)
        self.expectNCBI.addItems(str(i) for i in self.expectL_NCBI)
        self.expectNCBI.setCurrentIndex(10)
        self.filterNCBI.addItems(self.yesNoL)
        self.filterNCBI.setCurrentIndex(1)
        self.recordsNCBI.addItems(str(i) for i in self.scoreL_NCBI)
        self.recordsNCBI.setCurrentIndex(3)
        self.gapOpenNCBI.addItems(str(i) for i in self.gapOpenBL62L)
        self.gapOpenNCBI.setCurrentIndex(5)
        self.gapExtNCBI.addItems(str(i) for i in self.gapExtL_NCBI)
        self.gapExtNCBI.setCurrentIndex(1)
        self.gapAlign.addItems(self.yesNoL)
        self.gapAlign.setCurrentIndex(0)
        self.dropNCBI.addItems(str(i) for i in self.DropoffL)
        self.dropNCBI.setCurrentIndex(0)
        self.alignmentsNCBI.addItems(str(i) for i in self.scoreL_NCBI)
        self.alignmentsNCBI.setCurrentIndex(3)
        self.composit.addItems(str(i) for i in self.compositionalL)
        self.composit.setCurrentIndex(0)
        # ___ FASTA parameter value ___
        self.matrixTypesFASTA.addItems(self.matrixL_FASTA)
        self.matrixTypesFASTA.currentTextChanged.connect(self.get_gap)
        self.matrixTypesFASTA.setCurrentIndex(0)
        self.expectFASTA.addItems(str(i) for i in self.expectL_FASTA)
        self.expectFASTA.setCurrentIndex(0)
        self.expectUp.addItems(str(i) for i in self.expectUpL)
        self.expectUp.setCurrentIndex(11)
        self.filterFASTA.addItems(self.filterL_FASTA)
        self.filterFASTA.setCurrentIndex(0)
        self.recordsFASTA.addItems(str(i) for i in self.scoreL_FASTA)
        self.recordsFASTA.setCurrentIndex(4)
        self.statEst.addItems(self.statisticL)
        self.statEst.setCurrentIndex(0)
        self.gapOpenFASTA.addItems(str(i) for i in self.gapOpenL_FASTA)
        self.gapOpenFASTA.setCurrentIndex(10)
        self.gapExtFASTA.addItems(str(i) for i in self.gapExtL_FASTA)
        self.gapExtFASTA.setCurrentIndex(2)
        self.ktup.addItems(str(i) for i in self.ktupL)
        self.ktup.setCurrentIndex(1)
        self.alignmentsFASTA.addItems(str(i) for i in self.scoreL_FASTA)
        self.alignmentsFASTA.setCurrentIndex(4)
        self.progTypes.addItems(self.programL)
        # ___ FASTM parameter value ___
        self.matrixTypesFASTM.addItems(self.matrixL_FASTM)
        self.matrixTypesFASTM.currentTextChanged.connect(self.get_gap)
        self.matrixTypesFASTM.setCurrentIndex(6)
        self.expectFASTM.addItems(str(i) for i in self.expectL_FASTA)
        self.expectFASTM.setCurrentIndex(0)
        self.expectUpFASTM.addItems(str(i) for i in self.expectUpL)
        self.expectUpFASTM.setCurrentIndex(11)
        self.filterFASTM.addItems(self.filterL_FASTA)
        self.filterFASTM.setCurrentIndex(0)
        self.recordsFASTM.addItems(str(i) for i in self.scoreL_FASTA)
        self.recordsFASTM.setCurrentIndex(4)
        self.statEstFASTM.addItems(self.statisticL)
        self.statEstFASTM.setCurrentIndex(0)
        self.gapOpenFASTM.addItems(str(i) for i in self.gapOpenL_FASTM)
        self.gapOpenFASTM.setCurrentIndex(22)
        self.gapExtFASTM.addItems(str(i) for i in self.gapExtL_FASTM)
        self.gapExtFASTM.setCurrentIndex(4)
        self.ktupFASTM.addItems(str(i) for i in self.ktupL)
        self.ktupFASTM.setCurrentIndex(1)
        self.alignmentsFASTM.addItems(str(i) for i in self.scoreL_FASTA)
        self.alignmentsFASTM.setCurrentIndex(4)
        self.progTypesFASTM.addItems(self.programL_fastm)
        # ___ psi blast parameter value ___
        self.pssmPSIBL.addItems(str(i) for i in self.pssmL)
        self.pssmPSIBL.setCurrentIndex(0)
        self.matrixTypesPSIBL.addItems(self.matrixL_PSIBL)
        self.matrixTypesPSIBL.setCurrentIndex(0)
        self.matrixTypesPSIBL.currentTextChanged.connect(self.get_gap)
        self.expectPSIBL.addItems(str(i) for i in self.expectL_PSIBL)
        self.expectPSIBL.setCurrentIndex(10)
        self.filterPSIBL.addItems(self.yesNoL)
        self.filterPSIBL.setCurrentIndex(1)
        self.recordsPSIBL.addItems(str(i) for i in self.scoreL_PSIBL)
        self.recordsPSIBL.setCurrentIndex(8)
        self.dropPSIBL.addItems(str(i) for i in self.DropoffL_PSIBL)
        self.dropPSIBL.setCurrentIndex(6)
        self.fdrop.addItems(str(i) for i in self.DropoffL_PSIBL)
        self.fdrop.setCurrentIndex(8)
        self.gapOpenPSIBL.addItems(str(i) for i in self.gapOpenBL62L_PSIBL)
        self.gapOpenPSIBL.setCurrentIndex(4)
        self.gapExtPSIBL.addItems(str(i) for i in self.gapExtL_PSIBL)
        self.gapExtPSIBL.setCurrentIndex(2)
        self.alignmentsPSIBL.addItems(str(i) for i in self.scoreL_PSIBL)
        self.alignmentsPSIBL.setCurrentIndex(8)
        # ___ psi search parameter value ___
        self.pssmPSISE.addItems(str(i) for i in self.pssmL)
        self.pssmPSISE.setCurrentIndex(0)
        self.matrixTypesPSISE.addItems(self.matrixL_PSISE)
        self.matrixTypesPSISE.setCurrentIndex(0)
        self.matrixTypesPSISE.currentTextChanged.connect(self.get_gap)
        self.expectPSISE.addItems(str(i) for i in self.expectL_PSISE)
        self.expectPSISE.setCurrentIndex(10)
        self.filterPSISE.addItems(self.filterL_PSISE)
        self.filterPSISE.setCurrentIndex(0)
        self.recordsPSISE.addItems(str(i) for i in self.scoreL_PSISE)
        self.recordsPSISE.setCurrentIndex(8)
        self.alignmentsPSISE.addItems(str(i) for i in self.scoreL_PSIBL)
        self.alignmentsPSISE.setCurrentIndex(8)
        self.gapOpenPSISE.addItems(str(i) for i in self.gapOpenBL62L_PSISE)
        self.gapOpenPSISE.setCurrentIndex(4)
        self.gapExtPSISE_3.addItems(str(i) for i in self.gapExtL_PSISE)
        self.gapExtPSISE_3.setCurrentIndex(2)
        # search button
        self.searchBut.released.connect(self.get_params_th)

    def enterEvent(self, QEvent):
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']

        while self.path != path:
            self.path = path

    def browse_seq(self):
        try:
            seq_file = PyQt5.QtWidgets.QFileDialog.getOpenFileName()[0]
            with open(seq_file, 'r') as f:
                seq = ''.join(line.strip('\n') for line in f.readlines() if not line.startswith('>'))
                self.seq.setText(seq)
        except:
            pass

    def search_type(self):
        if self.searchType.currentText() == 'Whole sequence':
            self.residueRange.setEnabled(False)
            self.expect.setEnabled(False)
            self.residuerangeLabel.setEnabled(False)
            self.expectLabel.setEnabled(False)
        if self.searchType.currentText() == 'Residue range':
            self.residueRange.setEnabled(True)
            self.residuerangeLabel.setEnabled(False)
            self.residuerangeLabel.setEnabled(True)
            self.expectLabel.setEnabled(False)
            self.expect.setEnabled(False)
        if self.searchType.currentText() == 'Functional domain':
            self.residueRange.setEnabled(False)
            self.expect.setEnabled(True)
            self.residuerangeLabel.setEnabled(False)
            self.expectLabel.setEnabled(True)

    def show_blast_param(self):
        if self.blastParam.currentText() == 'NCBI BLASTp Parameters':
            self.group_NCBI.show()
            self.horizontalLayoutWidget_39.setGeometry(QtCore.QRect(20, 230, 951, 61))
            self.group_FASTA.hide()
            self.horizontalLayoutWidget_40.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_FASTM.hide()
            self.horizontalLayoutWidget_92.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_PSI_BLAST.hide()
            self.horizontalLayoutWidget_69.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_PSI_search.hide()
            self.horizontalLayoutWidget_88.setGeometry(QtCore.QRect(0, 0, 0, 0))
        elif self.blastParam.currentText() == 'FASTA BLAST Parameters':
            self.group_NCBI.hide()
            self.horizontalLayoutWidget_39.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_FASTA.show()
            self.horizontalLayoutWidget_40.setGeometry(QtCore.QRect(20, 230, 951, 61))
            self.group_FASTM.hide()
            self.horizontalLayoutWidget_92.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_PSI_BLAST.hide()
            self.horizontalLayoutWidget_69.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_PSI_search.hide()
            self.horizontalLayoutWidget_88.setGeometry(QtCore.QRect(0, 0, 0, 0))
        elif self.blastParam.currentText() == 'FASTM BLAST Parameters':
            self.group_NCBI.hide()
            self.horizontalLayoutWidget_39.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_FASTA.hide()
            self.horizontalLayoutWidget_40.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_FASTM.show()
            self.horizontalLayoutWidget_92.setGeometry(QtCore.QRect(20, 230, 951, 61))
            self.group_PSI_BLAST.hide()
            self.horizontalLayoutWidget_69.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_PSI_search.hide()
            self.horizontalLayoutWidget_88.setGeometry(QtCore.QRect(0, 0, 0, 0))
        elif self.blastParam.currentText() == 'PSI-BLAST Parameters':
            self.group_NCBI.hide()
            self.horizontalLayoutWidget_39.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_FASTA.hide()
            self.horizontalLayoutWidget_40.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_FASTM.hide()
            self.horizontalLayoutWidget_92.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_PSI_BLAST.show()
            self.horizontalLayoutWidget_69.setGeometry(QtCore.QRect(20, 230, 951, 61))
            self.group_PSI_search.hide()
            self.horizontalLayoutWidget_88.setGeometry(QtCore.QRect(0, 0, 0, 0))
        elif self.blastParam.currentText() == 'PSI-Search Parameters':
            self.group_NCBI.hide()
            self.horizontalLayoutWidget_39.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_FASTA.hide()
            self.horizontalLayoutWidget_40.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_FASTM.hide()
            self.horizontalLayoutWidget_92.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_PSI_BLAST.hide()
            self.horizontalLayoutWidget_69.setGeometry(QtCore.QRect(0, 0, 0, 0))
            self.group_PSI_search.show()
            self.horizontalLayoutWidget_88.setGeometry(QtCore.QRect(20, 230, 951, 61))

    def get_params_th(self):
        t = threading.Thread(target=self.get_params)
        t.start()
        self.searchBut.setEnabled(False)

    def get_params(self):
        self.msgLabel.setStyleSheet('color: #000000')
        self.msgLabel.setText('processing . . .')
        QApplication.processEvents()
        self.accessionL = []
        self.blastSelectedL = []
        self.search_result = {}
        self.domainDic = {}
        # title
        self.params['title'] = self.title.text()
        # expect
        self.params['evalue'] = float(self.expect.currentText())
        # program type
        for i in self.blastType.selectedItems():
            self.blastSelectedL.append(i.text())
        # database
        self.params['database'] = 'Pfam-A'.lower()
        # active site predict
        self.params['asp'] = 'false'
        # format
        self.params['format'] = 'json'
        # email
        self.params['email'] = 'pymodel.v1@gmail.com'
        # accession
        r = re.compile(r'(\S{4,})[\s+]?')
        s = r.findall(self.accession.text().strip())
        for i in s:
            if i not in self.accessionL:
                self.accessionL.append(i.upper())
        # residue range
        if self.searchType.currentText() == 'Residue range':
            if self.residueRange.text().strip() == '':
                self.searchType.setCurrentText('Whole sequence')
            elif self.residueRange.text().strip() != '':
                r2 = re.compile(r'(\d+)-(\d+)')
                s2 = r2.findall(self.residueRange.text())
                for i in s2:
                    if i not in self.resRangeL:
                        self.resRangeL.append(i)
                if len(self.resRangeL) == 0:
                    self.searchType.setCurrentText('Whole sequence')
        self.get_sequence()
        self.run_blast()

    def get_params_ncbi(self):
        # email
        self.params_ncbi['email'] = 'pymodel.v1@gmail.com'
        # program
        self.params_ncbi['program'] = 'blastp'
        # seq type
        self.params_ncbi['stype'] = 'protein'
        # db
        self.params_ncbi['database'] = 'pdb'
        # matrix
        self.params_ncbi['matrix'] = self.matrixTypesNCBI.currentText()
        # expect
        self.params_ncbi['exp'] = self.expectNCBI.currentText()
        # filter
        if self.filterNCBI.currentText() == 'Yes':
            self.params['filter'] = 'T'
        else:
            self.params['filter'] = 'F'
        # record
        self.params_ncbi['scores'] = int(self.recordsNCBI.currentText())
        # gap open and extension => from get_gap function
        # compositional adjustments
        self.params_ncbi['compstats'] = str(self.composit.currentText())
        # dropoff
        self.params_ncbi['dropoff'] = int(self.dropNCBI.currentText())
        # gap align
        if self.gapAlign.currentText() == 'yes':
            self.params_ncbi['gapalign'] = 'true'
        else:
            self.params_ncbi['gapalign'] = 'false'
        # alignment
        self.params_ncbi['alignments'] = int(self.alignmentsNCBI.currentText())
        return self.params_ncbi

    def get_params_fasta(self):
        # email
        self.params_fasta['email'] = 'pymodel.v1@gmail.com'
        # program
        self.params_fasta['program'] = self.progTypes.currentText().lower()
        # seq type
        self.params_fasta['stype'] = 'protein'
        # db
        self.params_fasta['database'] = 'pdb'
        # matrix
        for i in self.matrixL_FASTA:
            if self.matrixTypesFASTA.currentText() == i:
                self.params_fasta['matrix'] = str(self.matrixAbbrL_FASTA[self.matrixL_FASTA.index(i)])
        # expect up
        self.params_fasta['expupperlim'] = float(self.expectUp.currentText())
        # expect
        self.params_fasta['explowlim'] = float(self.expect.currentText())
        # filter
        self.params_fasta['filter'] = str(self.filterFASTA.currentText()).lower()
        # record
        self.params_fasta['scores'] = int(self.recordsFASTA.currentText())
        # gap open and extension => from get_gap function
        # annotation feature
        self.params_fasta['annotfeats'] = 'true'
        # statistical est
        for i in self.statisticL:
            if self.statEst.currentText() == i:
                self.params_fasta['stats'] = str(self.statisticAbbrL[self.statisticL.index(i)])
        # ktup
        for i in self.ktupL:
            if self.ktup.currentText() == i:
                self.params_fasta['ktup'] = int(self.ktupAbbrL[self.ktupL.index(i)])
        # hsp
        self.params_fasta['hsps'] = 'false'
        # histogram
        self.params_fasta['hist'] = 'false'
        # alignment
        self.params_fasta['alignments'] = int(self.alignmentsFASTA.currentText())
        # format
        self.params_fasta['scoreformat'] = '9'
        return self.params_fasta

    def get_params_fastm(self):
        # email
        self.params_fastm['email'] = 'pymodel.v1@gmail.com'
        # program
        self.params_fastm['program'] = self.progTypesFASTM.currentText().lower()
        # seq type
        self.params_fastm['stype'] = 'protein'
        # db
        self.params_fastm['database'] = 'pdb'
        # matrix
        for i in self.matrixL_FASTM:
            if self.matrixTypesFASTM.currentText() == i:
                self.params_fastm['matrix'] = str(self.matrixAbbrL_FASTM[self.matrixL_FASTM.index(i)])
        # expect up
        self.params_fastm['expupperlim'] = float(self.expectUp.currentText())
        # expect
        self.params_fastm['explowlim'] = float(self.expect.currentText())
        # filter
        self.params_fastm['filter'] = str(self.filterFASTM.currentText()).lower()
        # record
        self.params_fastm['scores'] = int(self.recordsFASTM.currentText())
        # gap open and extension => from get_gap function
        # annotation feature
        self.params_fastm['annotfeats'] = 'true'
        # statistical est
        for i in self.statisticL:
            if self.statEst.currentText() == i:
                self.params_fastm['stats'] = str(self.statisticAbbrL[self.statisticL.index(i)])
        # ktup
        for i in self.ktupL:
            if self.ktup.currentText() == i:
                self.params_fastm['ktup'] = int(self.ktupAbbrL[self.ktupL.index(i)])
        # hsp
        self.params_fastm['hsps'] = 'false'
        # histogram
        self.params_fastm['hist'] = 'false'
        # alignment
        self.params_fastm['alignments'] = int(self.alignmentsFASTM.currentText())
        # format
        self.params_fastm['scoreformat'] = '9'
        return self.params_fastm

    def get_params_psibl(self):
        # email
        self.params_psi_blast['email'] = 'pymodel.v1@gmail.com'
        # db
        self.params_psi_blast['database'] = 'pdb'
        # PSSM
        self.params_psi_blast['psithr'] = float(self.pssmPSIBL.currentText())
        # matrix
        self.params_psi_blast['matrix'] = self.matrixTypesPSIBL.currentText()
        # expect
        self.params_psi_blast['expthr'] = int(self.expect.currentText())
        # filter
        if self.filterPSIBL.currentText() == 'Yes':
            self.params_psi_blast['filter'] = 'T'
        else:
            self.params_psi_blast['filter'] = 'F'
        # record
        self.params_psi_blast['scores'] = int(self.recordsPSIBL.currentText())
        # gap open and extension => from get_gap function
        # dropoff
        self.params_psi_blast['dropoff'] = int(self.dropPSIBL.currentText())
        # final dropoff
        self.params_psi_blast['finaldropoff'] = int(self.fdrop.currentText())
        # alignment
        self.params_psi_blast['alignments'] = int(self.alignmentsPSIBL.currentText())
        return self.params_psi_blast

    def get_params_psise(self):
        # email
        self.params_psi_search['email'] = 'pymodel.v1@gmail.com'
        # db
        self.params_psi_search['database'] = 'pdb'
        # PSSM
        self.params_psi_search['psithr'] = float(self.pssmPSISE.currentText())
        # matrix
        self.params_psi_search['matrix'] = self.matrixTypesPSISE.currentText()
        # expect
        self.params_psi_search['expthr'] = int(self.expect.currentText())
        # filter
        self.params_psi_search['filter'] = self.filterPSISE.currentText().lower()
        # record
        self.params_psi_search['scores'] = int(self.recordsPSISE.currentText())
        # gap open and extension => from get_gap function
        # score format
        self.params_psi_search['scoreformat'] = 'default'
        # HOE mask
        self.params_psi_search['mask'] = 'true'
        # annotation Features
        self.params_psi_search['annotfeats'] = 'true'
        # alignment
        self.params_psi_search['alignments'] = int(self.alignmentsPSISE.currentText())
        return self.params_psi_search

    def get_gap(self):
        # ncbi
        if self.matrixTypesNCBI.currentText() == 'BLOSUM62':
            self.gapOpenNCBI.setCurrentIndex(5)
            self.gapExtNCBI.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenNCBI.currentText()
            self.params['gapext'] = self.gapExtNCBI.currentText()
        elif self.matrixTypesNCBI.currentText() == 'PAM30':
            self.gapOpenNCBI.setCurrentIndex(4)
            self.gapExtNCBI.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenNCBI.currentText()
            self.params['gapext'] = self.gapExtNCBI.currentText()
        elif self.matrixTypesNCBI.currentText() == 'PAM70':
            self.gapOpenNCBI.setCurrentIndex(4)
            self.gapExtNCBI.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenNCBI.currentText()
            self.params['gapext'] = self.gapExtNCBI.currentText()
        elif self.matrixTypesNCBI.currentText() == 'PAM250':
            self.gapOpenNCBI.setCurrentIndex(3)
            self.gapExtNCBI.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenNCBI.currentText()
            self.params['gapext'] = self.gapExtNCBI.currentText()
        elif self.matrixTypesNCBI.currentText() == 'BLOSUM80':
            self.gapOpenNCBI.setCurrentIndex(4)
            self.gapExtNCBI.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenNCBI.currentText()
            self.params['gapext'] = self.gapExtNCBI.currentText()
        elif self.matrixTypesNCBI.currentText() == 'BLOSUM45':
            self.gapOpenNCBI.setCurrentIndex(4)
            self.gapExtNCBI.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenNCBI.currentText()
            self.params['gapext'] = self.gapExtNCBI.currentText()
        elif self.matrixTypesNCBI.currentText() == 'BLOSUM50':
            self.gapOpenNCBI.setCurrentIndex(4)
            self.gapExtNCBI.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenNCBI.currentText()
            self.params['gapext'] = self.gapExtNCBI.currentText()
        elif self.matrixTypesNCBI.currentText() == 'BLOSUM90':
            self.gapOpenNCBI.setCurrentIndex(5)
            self.gapExtNCBI.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenNCBI.currentText()
            self.params['gapext'] = self.gapExtNCBI.currentText()
        # fasta
        if self.matrixTypesFASTA.currentText() == 'BLOSUM50':
            self.gapOpenFASTA.setCurrentIndex(10)
            self.gapExtFASTA.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'BLASTP62':
            self.gapOpenFASTA.setCurrentIndex(11)
            self.gapExtFASTA.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'BLOSUM80':
            self.gapOpenFASTA.setCurrentIndex(10)
            self.gapExtFASTA.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'PAM250':
            self.gapOpenFASTA.setCurrentIndex(10)
            self.gapExtFASTA.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'PAM120':
            self.gapOpenFASTA.setCurrentIndex(16)
            self.gapExtFASTA.setCurrentIndex(4)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'MDM40':
            self.gapOpenFASTA.setCurrentIndex(12)
            self.gapExtFASTA.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'MDM20':
            self.gapOpenFASTA.setCurrentIndex(22)
            self.gapExtFASTA.setCurrentIndex(4)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'MDM10':
            self.gapOpenFASTA.setCurrentIndex(23)
            self.gapExtFASTA.setCurrentIndex(4)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'VTML160':
            self.gapOpenFASTA.setCurrentIndex(12)
            self.gapExtFASTA.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'VTML120':
            self.gapOpenFASTA.setCurrentIndex(11)
            self.gapExtFASTA.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'VTML80':
            self.gapOpenFASTA.setCurrentIndex(11)
            self.gapExtFASTA.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'VTML40':
            self.gapOpenFASTA.setCurrentIndex(12)
            self.gapExtFASTA.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'VTML20':
            self.gapOpenFASTA.setCurrentIndex(15)
            self.gapExtFASTA.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        elif self.matrixTypesFASTA.currentText() == 'VTML10':
            self.gapOpenFASTA.setCurrentIndex(16)
            self.gapExtFASTA.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTA.currentText()
            self.params['gapext'] = self.gapExtFASTA.currentText()
        # program - ktup check
        # if self.progTypes.currentText() == 'FASTA':
        #     self.ktup.setCurrentIndex(4)
        # else:
        #     self.ktup.setCurrentIndex(6)
        # fastm
        if self.matrixTypesFASTM.currentText() == 'BLOSUM50':
            self.gapOpenFASTM.setCurrentIndex(10)
            self.gapExtFASTM.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
        elif self.matrixTypesFASTM.currentText() == 'BLOSUM62':
            self.gapOpenFASTM.setCurrentIndex(7)
            self.gapExtFASTM.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
        elif self.matrixTypesFASTM.currentText() == 'BLASTP62':
            self.gapOpenFASTM.setCurrentIndex(11)
            self.gapExtFASTM.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
        elif self.matrixTypesFASTM.currentText() == 'BLOSUM80':
            self.gapOpenFASTM.setCurrentIndex(16)
            self.gapExtFASTM.setCurrentIndex(4)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
        elif self.matrixTypesFASTM.currentText() == 'PAM250':
            self.gapOpenFASTM.setCurrentIndex(10)
            self.gapExtFASTM.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
        elif self.matrixTypesFASTM.currentText() == 'PAM120':
            self.gapOpenFASTM.setCurrentIndex(16)
            self.gapExtFASTM.setCurrentIndex(4)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
        elif self.matrixTypesFASTM.currentText() == 'MDM40':
            self.gapOpenFASTM.setCurrentIndex(12)
            self.gapExtFASTM.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
        elif self.matrixTypesFASTM.currentText() == 'MDM20':
            self.gapOpenFASTM.setCurrentIndex(22)
            self.gapExtFASTM.setCurrentIndex(4)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
        elif self.matrixTypesFASTM.currentText() == 'MDM10':
            self.gapOpenFASTM.setCurrentIndex(23)
            self.gapExtFASTM.setCurrentIndex(4)
            self.params['gapopen'] = self.gapOpenFASTM.currentText()
            self.params['gapext'] = self.gapExtFASTM.currentText()
            # psi blast
        if self.matrixTypesPSIBL.currentText() == 'BLOSUM45':
            self.gapOpenPSIBL.setCurrentIndex(4)
            self.gapExtPSIBL.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenPSIBL.currentText()
            self.params['gapext'] = self.gapExtPSIBL.currentText()
        elif self.matrixTypesPSIBL.currentText() == 'BLOSUM62':
            self.gapOpenPSIBL.setCurrentIndex(0)
            self.gapExtPSIBL.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenPSIBL.currentText()
            self.params['gapext'] = self.gapExtPSIBL.currentText()
        elif self.matrixTypesPSIBL.currentText() == 'BLOSUM80':
            self.gapOpenPSIBL.setCurrentIndex(2)
            self.gapExtPSIBL.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenPSIBL.currentText()
            self.params['gapext'] = self.gapExtPSIBL.currentText()
        elif self.matrixTypesPSIBL.currentText() == 'PAM30':
            self.gapOpenPSIBL.setCurrentIndex(1)
            self.gapExtPSIBL.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenPSIBL.currentText()
            self.params['gapext'] = self.gapExtPSIBL.currentText()
        elif self.matrixTypesPSIBL.currentText() == 'PAM70':
            self.gapOpenPSIBL.setCurrentIndex(2)
            self.gapExtPSIBL.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenPSIBL.currentText()
            self.params['gapext'] = self.gapExtPSIBL.currentText()
        # psi search
        if self.matrixTypesPSISE.currentText() == 'BLOSUM45':
            self.gapOpenPSISE.setCurrentIndex(4)
            self.gapExtPSISE_3.setCurrentIndex(2)
            self.params['gapopen'] = self.gapOpenPSISE.currentText()
            self.params['gapext'] = self.gapExtPSISE_3.currentText()
        elif self.matrixTypesPSISE.currentText() == 'BLOSUM62':
            self.gapOpenPSISE.setCurrentIndex(3)
            self.gapExtPSISE_3.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenPSISE.currentText()
            self.params['gapext'] = self.gapExtPSISE_3.currentText()
        elif self.matrixTypesPSISE.currentText() == 'BLOSUM80':
            self.gapOpenPSISE.setCurrentIndex(2)
            self.gapExtPSISE_3.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenPSISE.currentText()
            self.params['gapext'] = self.gapExtPSISE_3.currentText()
        elif self.matrixTypesPSISE.currentText() == 'PAM30':
            self.gapOpenPSISE.setCurrentIndex(1)
            self.gapExtPSISE_3.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenPSISE.currentText()
            self.params['gapext'] = self.gapExtPSISE_3.currentText()
        elif self.matrixTypesPSISE.currentText() == 'PAM70':
            self.gapOpenPSISE.setCurrentIndex(2)
            self.gapExtPSISE_3.setCurrentIndex(1)
            self.params['gapopen'] = self.gapOpenPSISE.currentText()
            self.params['gapext'] = self.gapExtPSISE_3.currentText()

    # _______ analysis _______

    def progressbar(self, value, endvalue, bar_length):
        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length) - 1) + '>'
        spaces = ' ' * (bar_length - len(arrow))
        sys.stdout.write("\r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

    def get_sequence(self):
        try:
            self.msgLabel.setStyleSheet('color: #000000')
            self.msgLabel.setText('Retrieving sequence . . .')
            QApplication.processEvents()
            if len(self.accessionL) > 0:
                for acc in self.accessionL:
                    self.get_seq_from_servers(acc)
            if len(self.seq.toPlainText().strip()) > 0:
                self.search_result['sequence'] = [self.seq.toPlainText().strip(), 'your_sequence']
        except requests.exceptions.ConnectionError:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Connection error')
            print('Connection error.')
        except Exception as er:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Error')
            self.searchBut.setEnabled(True)
            print(er)

    def get_seq_from_servers(self, acc):
        # get path
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']
        # connect to uniprot
        if len(acc) >= 6:
            try:
                url = 'https://www.uniprot.org/uniprot/' + acc + '.xml'
                req = requests.get(url)
                req.raise_for_status()
                soup = bs4.BeautifulSoup(req.content, 'lxml')
                accSwTr = []
                entrySw = soup.find_all('entry', {'dataset': 'Swiss-Prot'})
                entryTr = soup.find_all('entry', {'dataset': 'TrEMBL'})
                for i in entrySw:
                    i['id'] = 'SwpTr'
                for i in entryTr:
                    i['id'] = 'SwpTr'
                entrySwTr = soup.find_all('entry', id='SwpTr', limit=1)
                for ac in entrySwTr:
                    accSwTr.append(ac.accession.getText())
                # accNum = accSwTr
                if len(accSwTr) > 0:
                    for n in range(len(accSwTr)):
                        sequence = entrySwTr[n].find_all('sequence')
                        fullnameSeqSwTr = entrySwTr[n].find('protein')
                        for seq in sequence:
                            self.search_result[acc] = [seq.getText().strip()]
                        self.search_result[accSwTr[n]].insert(1, fullnameSeqSwTr.fullname.getText().strip())
            except:
                # entrez search
                Entrez.email = self.params['email']
                handle = Entrez.efetch(db='protein', id=acc, rettype='gb')
                handleRead = SeqIO.read(handle, 'gb')
                self.search_result[acc] = [handleRead.seq.__str__(), handleRead.name]
        elif len(acc) == 4:
            # get pdb file
            pdb = PDBList()
            pdbfile = pdb.retrieve_pdb_file(acc, file_format='pdb', obsolete=False, pdir=path)
            # get name
            p = PDBParser()
            s = p.get_structure(acc, os.path.join(path, acc + '.pdb'))
            # pdb to fasta
            st = set()
            readPDB = open(pdbfile, "rU")
            for record in SeqIO.parse(readPDB, "pdb-seqres"):
                st.add(record.seq)
            self.search_result[acc] = [''.join([str(i) for i in st]), s.header['name']]
            readPDB.close()
            os.unlink(os.path.join(path, acc.lower() + '.pdb'))

    def run_blast(self):
        try:
            cycle = len(self.blastSelectedL)
            for k, v in self.search_result.items():
                if self.searchType.currentText() == 'Functional domain':
                    self.archive()
                    self.msgLabel.setStyleSheet('color: #000000')
                    self.msgLabel.setText('Domain analysis . . .')
                    QApplication.processEvents()
                    print('Domain analysis for', str('\"' + k + '\"') + '...')
                    self.params['sequence'] = v[0]
                    run_pfamscan.runPfamscan.runpfamscan(self.params)
                    run_pfamscan.runPfamscan.getResult(k)
                    self.domainanalysis(v[0])
                    print('Finished.\n')
                elif self.searchType.currentText() == 'Residue range':
                    self.archive()
                    self.get_ranges(seq=v[0], acc=k)
                elif self.searchType.currentText() == 'Whole sequence':
                    self.archive()
                    self.domainDic[k] = [1, len(v[0]), v[0]]
                if len(self.domainDic) > 0:
                    # create target pir file
                    self.target_pir_file(k, v)
                    # run blast
                    for program in self.blastSelectedL:
                        self.msgLabel.setStyleSheet('color: #000000')
                        self.msgLabel.setText('Running ' + program + ' . . .')
                        QApplication.processEvents()
                        self.blast(acc=k, program=program)
                        if cycle > 1:
                            self.archive()
                        cycle -= 1
                        # domainAna.domainDic.clear()
                elif len(self.domainDic) == 0:
                    self.msgLabel.setStyleSheet('color: red')
                    self.msgLabel.setText('Sequence not found.')
                    QApplication.processEvents()
                    print('Sequence not found.')
                cycle = len(self.blastSelectedL)
            if len(self.domainDic) > 0:
                self.msgLabel.setStyleSheet('color: green')
                self.msgLabel.setText('Finished')
                self.searchBut.setEnabled(True)
            else:
                self.msgLabel.setText('')
        except requests.exceptions.ConnectionError:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Connection error')
            print('Connection error.')
        except Exception as er:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Error')
            self.searchBut.setEnabled(True)
            print(er)

    def archive(self):
        # get accession as name of archive dir
        r = re.compile(r'Accession: (.*)')
        s = ' '
        files = []
        num = 1
        try:
            for i in os.listdir(self.path):
                if i.endswith('.txt'):
                    file = open(os.path.join(self.path, i), 'r').read()
                    s = r.search(file)
                    files = [i for i in os.listdir(self.path) if os.path.isfile(os.path.join(self.path, i))]
            # create archive dir
            if len(files) > 0:
                while True:
                    dirName = s.group(1) + '_' + str(num)
                    if not os.path.exists(os.path.join(self.path, dirName)):
                        os.mkdir(os.path.join(self.path, dirName))
                    else:
                        num += 1
                        continue
                    for i in os.listdir(self.path):
                        if os.path.isfile(os.path.join(self.path, i)):
                            shutil.move(os.path.join(self.path, i), os.path.join(self.path, dirName))
                    break
        except:
            pass

    def domainanalysis(self, seq):
        # make domain name as keys & range and seq as values of dic
        domCount = 0
        num = 1
        frm = int()
        to = int()
        for file in os.listdir(self.path):
            if file.endswith('.domain'):
                fileRead = open(os.path.join(self.path, file), 'r')
                jsonLoad = json.load(fileRead)
                for n in range(len(jsonLoad)):
                    for k, v in jsonLoad[n].items():
                        if k == 'env':
                            domCount += 1
                            frm = v['from']
                            to = v['to']
                        if k == 'desc':
                            self.domainDic[str(num) + ') ' + v] = [int(frm), int(to)]
                    num += 1
                    fileRead.close()
        # add seq
        for k, v in self.domainDic.items():
            v.insert(2, seq[v[0] - 1: v[1]])
        # check to get all domain
        while domCount > len(self.domainDic):
            self.domainanalysis(seq=seq)

    def get_ranges(self, seq, acc):
        for i in self.resRangeL:
            self.domainDic['seq_' + '(' + str(i[0]) + '-' + str(i[1]) + ')' + '(' + acc + ')'] = [int(i[0]), int(i[1]),
                            seq[int(i[0]) - 1: int(i[1])]]
            with open(os.path.join(self.path, acc + '.domain'), 'w') as f:
                for k, v in self.domainDic.items():
                    f.write(k + ': ' + v[2] + '\n\n')

    def target_pir_file(self, name, seq):
        # create pir format file
        pirFile = name + '.pir'
        pirCode = 'P1'
        with open(os.path.join(self.path, pirFile), 'w') as f:
            f.write(
                '>' + pirCode + ';' + name + '\n' + 'sequence:' + name + ':::::::0.00: 0.00' + '\n' + seq[0] + '*')
        # add file name to json file
        with open('info', 'r+') as f:
            jFile = json.load(f)
            jFile['Target'] = name
            f.seek(0)
            f.write(json.dumps(jFile))
            f.truncate()

    def blast(self, acc, program):
        # create info file
        infoName = acc + '_' + program + '_info' + '.txt'
        with open(os.path.join(self.path, infoName), 'a') as f:
            f.write(
                'Title: ' + self.params[
                    'title'] + '\n' + 'Accession: ' + acc + '\n' + 'Program: ' + program + '\n' + 'Job IDs: \n')

        # ___run ncbi blast___
        if program == 'NCBI BLASTp':
            print('Running NCBI BLASTp...')
            num = 1
            for k, v in self.domainDic.items():
                self.progressbar(int(num), int(len(self.domainDic)), 20)
                self.params_ncbi['sequence'] = v[2]
                run_ncbi_blast.runNcbiBlast.runncbiblast(self.get_params_ncbi())
                run_ncbi_blast.runNcbiBlast.getResult(k)
                # write job ids
                with open(os.path.join(self.path, infoName), 'a') as f:
                    f.write('-  ' + k + ' (' + str(v[0]) + ':' + str(v[1]) + ')' + ': ' + run_ncbi_blast.runNcbiBlast.jobId)
                    f.write('\n\n')
                num += 1
                time.sleep(5)
            print('\nFinished\n\n')

        # ___run fasta blast___
        elif program == 'FASTA BLAST':
            print('Running FASTA BLAST...')
            num = 1
            for k, v in self.domainDic.items():
                self.progressbar(int(num), int(len(self.domainDic)), 20)
                self.params_fasta['sequence'] = v[2]
                run_fasta_blast.runFastaBlast.runfasta(self.get_params_fasta())
                run_fasta_blast.runFastaBlast.getResult(k)
                # write job ids
                with open(os.path.join(self.path, infoName), 'a') as f:
                    f.write('-  ' + k + ' (' + str(v[0]) + ':' + str(v[1]) + ')' + ': ' + run_fasta_blast.runFastaBlast.jobId)
                    f.write('\n\n')
                num += 1
                time.sleep(5)
            print('\nFinished\n\n')

        # ___run psi-blast___
        elif program == 'PSI-BLAST':
            print('Running PSI-BLAST...')
            num = 1
            for k, v in self.domainDic.items():
                self.progressbar(int(num), int(len(self.domainDic)), 20)
                self.params_psi_blast['sequence'] = v[2]
                run_psi_blast.runPsiBlast.runpsiblast(self.get_params_psibl())
                run_psi_blast.runPsiBlast.getResult(k)
                # write job ids
                with open(os.path.join(self.path, infoName), 'a') as f:
                    f.write('-  ' + k + ' (' + str(v[0]) + ':' + str(v[1]) + ')' + ': ' + run_psi_blast.runPsiBlast.jobId)
                    f.write('\n\n')
                num += 1
                time.sleep(5)
            print('\nFinished\n\n')

            # ___run psi-search___
        elif program == 'PSI-Search':
            print('Running PSI-Search...')
            num = 1
            for k, v in self.domainDic.items():
                self.progressbar(int(num), int(len(self.domainDic)), 20)
                self.params_psi_search['sequence'] = v[2]
                run_psi_search.runPsisearch.runpsisearch(self.get_params_psise())
                run_psi_search.runPsisearch.getResult(k)
                # write job ids
                with open(os.path.join(self.path, infoName), 'a') as f:
                    f.write('-  ' + k + ' (' + str(v[0]) + ':' + str(v[1]) + ')' + ': ' + run_psi_search.runPsisearch.jobId)
                    f.write('\n\n')
                f.close()
                num += 1
                time.sleep(5)
            print('\nFinished\n\n')


# def main():
#     app = QApplication(sys.argv)
#     form = TempSearch()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()