import os
import re
import sys
import json
import shutil
import threading
import bs4
import requests
from Bio import Entrez
from Bio import SeqIO
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog
import tpl_align_temp
import view_align
import config


class AlignTemp(QMainWindow, tpl_align_temp.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    pir_errorL = []
    seqL = []
    pdb_fileL = []
    pir_fileL = []
    ali_file = ''
    YesNoL = ['Yes', 'No']
    patchTypeL = [['1MC', '+'], ['1MT', '+'], ['25P1', '+'], ['25P2', '+'], ['3PHO', '+'], ['3TER', '+'], ['5DP', '+'],
                ['5MC1', '+'], ['5MC2', '+'], ['5MET', '+'], ['5PHO', '+'], ['5TER', '+'], ['9EG', '+'],
                ['9MA', '+'], ['9MG', '+'], ['ACE', 1], ['ACP', 1], ['ASPP', '+'], ['CT1', 1], ['CT2', 1],
                ['CT3', 1], ['CTER', 1], ['DEO1', '+'], ['DEO2', '+'], ['DISU', 2], ['FHEM', '+'], ['GLUP', '+'],
                ['GLYP', 1], ['INO1', '+'], ['LIG1', 2], ['LIG2', 2], ['LIG3', 2], ['LINK', 2], ['NTER', 1],
                ['PHEM', 2], ['PLIG', 3], ['PLO2', 3], ['PROP', 1], ['PURA', '+'], ['PURG', '+'],
                ['PYRC', '+'], ['PYRU', '+'], ['TP1', '+'], ['TP1A', '+'], ['TP2', '+'], ['TP2A', '+']]
    # ___ alignment ___
    alignchainL = []
    alignaddedPdbChainDic = {}
    alignTypeL = ['Pairwise', 'Tree', 'Progressive']
    matrixNameL = ['BLOSUM62', 'AS1', 'PAM120', 'GONNET']
    matrixFileL = ['$(LIB)/blosum62.sim.mat', '$(LIB)/as1.sim.mat', '$(LIB)/120pam.sim.mat', '$(LIB)/gonnet.sim.mat']
    gap2dL = ['Helix', 'Beta', 'Accessibility', 'Straightness', 'CA-CA distance factor', 'DST min', 'DST power', 'T', 'Structure_profile']
    outputL = ['Alignment', 'Quality', 'Both']
    userMatrixTypeL = ['Similar', 'Distance']
    user_matrix_file = ''
    getLocalAlignL = []
    getMatrixL = []
    getGapFunctionL = []
    getAutoOverhangL = []
    getSimilarityFlagL = []
    getOutputL = []
    getDendrogramL = []
    getUsermatrixL = []
    featureWeightDic = {}
    getUseFitL = []
    getWriteFitL = []

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        #
        self.browse_ali.released.connect(self.browse_ali_file)
        #
        self.gen_ali.released.connect(self.get_seq_th)
        # ____ alignment ____
        # browse pdb files
        self.browsePDB.released.connect(self.browse_pdb_th)
        self.alignPDB.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        self.alignPDB.clicked.connect(self.align_get_pdb_chain)
        # add pdb/chain
        self.alignAddBut.released.connect(self.align_add_pdb_chain_th)
        # add all pdb/chain
        self.addAll.released.connect(self.add_all_pdb_th)
        # min pdb/chain
        self.alignMinBut.released.connect(self.align_min_pdb_chain)
        # align type
        self.alignType.addItems(self.alignTypeL)
        self.alignType.setCurrentIndex(1)
        # local alignment
        self.localGlobal.addItems(self.YesNoL)
        self.localGlobal.setCurrentIndex(1)
        self.localGlobal.currentIndexChanged.connect(self.check_local)
        # matrix
        self.matrix.addItems(self.matrixNameL)
        # gap function
        self.gapFunction.addItems(self.YesNoL)
        self.gapFunction.setCurrentIndex(1)
        # 1d gap open
        self.gapOpen1d.setText('-450')
        # 1d gap extension
        self.gapExten1d.setText('-50')
        # 2d helix gap penalty
        self.gapPenalt2dHelix.setText('3.5')
        # 2d beta gap penalty
        self.gapPenalt2dBeta.setText('3.5')
        # 2d access gap penalty
        self.gapPenalt2dAccess.setText('3.5')
        # 2d stright gap penalty
        self.gapPenalt2dStright.setText('0.2')
        # 2d ca-ca gap penalty
        self.gapPenalt2dCACA.setText('4.0')
        # 2d dst min gap penalty
        self.gapPenalt2dDstMin.setText('6.5')
        # 2d dat power gap penalty
        self.gapPenalt2dDstPower.setText('2.0')
        # 2d structure profile gap penalty
        self.gapPenalt2dStrucProfile.setText('0')
        # 3d gap open
        self.gapOpen3d.setText('0')
        # 3d gap extension
        self.gapExten3d.setText('3')
        # max gap length
        self.maxgapLen.setText('20')
        # gap-residue score
        self.gapResScore.setText('0')
        # gap-gap score
        self.gapGapScore.setText('0')
        # overhang
        self.overhang.setText('0')
        # auto-overhang
        self.autoOverhang.addItems(self.YesNoL)
        self.autoOverhang.setCurrentIndex(1)
        self.autoOverhang.currentIndexChanged.connect(self.check_auto_overhang)
        # overhang limit
        self.autoOverhangLim.setText('60')
        # overhang factor
        self.overhangFactor.setText('0.4')
        # similarity flag
        self.similarityFlag.addItems(self.YesNoL)
        self.similarityFlag.setCurrentIndex(1)
        # rms cut-off
        self.rmsCutOff.setText('3.5')
        # matrix offset
        self.matrixOffset.setText('0')
        # check user matrix
        self.userMatrix.released.connect(self.user_matrix)
        # matrix type
        self.userMatrixType.addItems(self.userMatrixTypeL)
        self.userMatrixType.setCurrentIndex(1)
        # use fit
        self.useFit.addItems(self.YesNoL)
        self.useFit.setCurrentIndex(1)
        # write fit
        self.writeFit.addItems(self.YesNoL)
        self.writeFit.setCurrentIndex(1)
        # dendrogram
        self.dendrogram.addItems(self.YesNoL)
        self.dendrogram.setCurrentIndex(1)
        # target sequnce
        self.seq.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.pir')])
        self.seq.setCurrentIndex(0)
        # check for pdb and ali files
        self.pdb_fileL = [os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')]
        self.pir_fileL = [os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pir')]
        # align button
        self.alignBut.released.connect(self.align_get_values_th)
        # view alignment
        self.viewAlign.released.connect(self.view_align)
        # remove file
        self.ali_remove_file.released.connect(self.remove_ali_file)

    def enterEvent(self, QEvent):
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']
        self.path = path
        #
        pdb_files = []
        pir_files = []
        for i in os.listdir(self.path):
            if i.endswith('.pdb') and i not in pdb_files:
                pdb_files.append(i)
            if i.endswith('.pir') and i not in pir_files:
                pir_files.append(i)
        while self.pdb_fileL != pdb_files:
            self.alignPDB.clear()
            self.alignPDB.addItems([os.path.basename(i) for i in os.listdir(path) if i.endswith('.pdb')])
            self.pdb_fileL = pdb_files
        while self.pir_fileL != pir_files:
            self.seq.clear()
            self.seq.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.pir')])
            self.pir_fileL = pir_files

    def browse_pdb_th(self):
        threading.Thread(target=self.browse_pdb).start()

    def browse_pdb(self):
        user_pdb = QFileDialog.getOpenFileNames()[0]
        if user_pdb:
            for file in user_pdb:
                shutil.copy(file, self.path)
            self.alignPDB.clear()
            self.alignPDB.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])

    def browse_ali_file(self):
        self.ali_file = QFileDialog.getOpenFileName()[0]
        self.browse_label.setText(os.path.basename(self.ali_file))

    def get_seq_th(self):
        t = threading.Thread(target=self.get_seq)
        t.start()
        self.gen_ali.setEnabled(False)

    def get_seq(self):
        try:
            self.pir_errorL = []
            self.seqL = []
            # entered seq
            enter_seq = self.seq_ali.toPlainText().strip()
            if len(enter_seq.strip()) > 0:
                if len(enter_seq.split('>')) == 1:
                    self.seqL.append(['your_seq', enter_seq.split('>')[0]])
                for n, s in enumerate(enter_seq.split('>')):
                    seqs = s.split('\n')
                    iden = re.findall(r'[\w\d]+', seqs[0])
                    if len(iden) == 0:
                        iden = 'your_seq'
                    elif len(iden) == 1:
                        iden = iden[0]
                    elif len(iden) > 1:
                        iden = iden[0] + '_' + iden[1]
                    seq = ''.join(i for i in seqs[1:])
                    if len(seq) > 0:
                        self.seqL.append([iden, seq])
            # entered acc
            acc = re.findall(r'(\S{4,})', self.acc_ali.text().strip())
            if len(acc) > 0:
                for n, i in enumerate(acc):
                    self.msg_gen.setStyleSheet('color: black')
                    self.msg_gen.setText('Retrieving %d/%d' % (n + 1, len(acc)))
                    QApplication.processEvents()
                    self.get_seq_from_servers(i.upper())
            # generate PIR files
            for name, seq in self.seqL:
                self.generate_ali_file(name, seq)
            if len(self.seqL) > 0 or len(acc) > 0:
                self.seq.clear()
                self.seq.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.pir')])
                self.msg_gen.setStyleSheet('color: green')
                self.msg_gen.setText('Finished')
            if len(self.pir_errorL) > 0:
                self.msg_gen.setStyleSheet('color: red')
                self.msg_gen.setText('Error in: ' + ', '.join(i for i in self.pir_errorL))
            QApplication.processEvents()
            self.gen_ali.setEnabled(True)
        except Exception as er:
            self.msg_gen.setStyleSheet('color: red')
            self.msg_gen.setText('Error')
            self.gen_ali.setEnabled(True)
            pass
            # print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))

    def get_seq_from_servers(self, acc):
        # connect to uniprot
        try:
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
                            for seq in sequence:
                                self.seqL.append([acc, seq.getText().strip()])
                except:
                    # entrez search
                    Entrez.email = 'pymodel.v1@gmail.com'
                    handle = Entrez.efetch(db='protein', id=acc, rettype='gb')
                    handleRead = SeqIO.read(handle, 'gb')
                    self.seqL.append([acc, handleRead.seq.__str__()])
            elif len(acc) == 4:
                # get pdb file
                pdb = PDBList()
                pdbfile = pdb.retrieve_pdb_file(acc, file_format='pdb', obsolete=False, pdir=self.path)
                os.rename(os.path.join(self.path, 'pdb' + acc.lower() + '.ent'), os.path.join(self.path, acc + '.pdb'))
                # get name
                p = PDBParser()
                s = p.get_structure(acc, os.path.join(self.path, acc + '.pdb'))
                # pdb to fasta
                st = set()
                readPDB = open(os.path.join(self.path, acc + '.pdb'), "rU")
                for record in SeqIO.parse(readPDB, "pdb-seqres"):
                    st.add(record.seq)
                self.seqL.append([acc, ''.join([str(i) for i in st])])
                readPDB.close()
                # os.unlink(os.path.join(self.path, acc.lower() + '.pdb'))
        except:
            self.pir_errorL.append(acc)

    def generate_ali_file(self, name, seq):
        # create pir format file
        pirFile = name + '.pir'
        pirCode = 'P1'
        with open(os.path.join(self.path, pirFile), 'w') as f:
            f.write('>' + pirCode + ';' + name + '\n' + 'sequence:' + name + ':::::::0.00: 0.00' + '\n' + seq + '*')

    def align_get_pdb_chain(self):
        try:
            self.alignchainL.clear()
            p = PDBParser()
            s = p.get_structure(self.alignPDB.currentItem().text(), os.path.join(self.path, self.alignPDB.currentItem().text()))
            for i in s.get_chains():
                self.alignchainL.append(i.id)
            self.alignChainPDB.clear()
            self.alignChainPDB.addItems(self.alignchainL)
            self.alignChainPDB.setCurrentRow(0)
        except:
            pass

    def align_add_pdb_chain_th(self):
        t = threading.Thread(target=self.align_add_pdb_chain)
        t.start()

    def align_add_pdb_chain(self):
        try:
            # add one file
            if len(self.alignPDB.selectedItems()) > 0 :
                self.alignAddBut.setText('Adding..')
                QApplication.processEvents()
                chainL = []
                self.alignpdbChainAdded.clear()
                selected = self.alignChainPDB.selectedItems()
                for i in selected:
                    chainL.append(i.text())

                k = [self.alignPDB.currentItem().text(), [i for i in chainL]]
                v0 = ''
                if len(chainL) == 0:
                    v0 = k[0]
                if len(chainL) == 1:
                    v0 = k[0] + ', ' + k[1][0]
                if len(chainL) > 1:
                    v0 = k[0] + ', ' + k[1][0] + ':' + k[1][len(k[1]) - 1]
                self.alignaddedPdbChainDic[k[0]] = [v0, os.path.join(self.path, k[0]), k[1]]
                self.alignpdbChainAdded.addItems([v[0] for v in self.alignaddedPdbChainDic.values()])
                # add multiple file
                p = PDBParser()
                if len(self.alignPDB.selectedItems()) > 1:
                    for file in self.alignPDB.selectedItems():
                        s = p.get_structure(file.text(), os.path.join(self.path, file.text()))
                        chainL = [i.id for i in s.get_chains() if i.id.strip() != '']

                        v0 = ''
                        if len(chainL) == 0:
                            v0 = file.text()
                        if len(chainL) == 1:
                            v0 = file.text() + ', ' + chainL[0]
                        if len(chainL) > 1:
                            v0 = file.text() + ', ' + chainL[0] + ':' + chainL[len(chainL) - 1]
                        self.alignaddedPdbChainDic[file.text()] = [v0, os.path.join(self.path, file.text()), chainL]
                        self.alignpdbChainAdded.clear()
                        self.alignpdbChainAdded.addItems([v[0] for v in self.alignaddedPdbChainDic.values()])
                self.alignAddBut.setText('Add')
                QApplication.processEvents()
        except:
            pass

    def add_all_pdb_th(self):
        t = threading.Thread(target=self.add_all_pdb)
        t.start()

    def add_all_pdb(self):
        try:
            self.addAll.setText('Adding..')
            QApplication.processEvents()
            #
            p = PDBParser()
            for file in [os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')]:
                s = p.get_structure(file, os.path.join(self.path, file))
                chainL = [i.id for i in s.get_chains() if i.id.strip() != '']

                v0 = ''
                if len(chainL) == 0:
                    v0 = file
                if len(chainL) == 1:
                    v0 = file + ', ' + chainL[0]
                if len(chainL) > 1:
                    v0 = file + ', ' + chainL[0] + ':' + chainL[len(chainL) - 1]
                self.alignaddedPdbChainDic[file] = [v0, os.path.join(self.path, file), chainL]
                self.alignpdbChainAdded.clear()
                self.alignpdbChainAdded.addItems([v[0] for v in self.alignaddedPdbChainDic.values()])
            self.addAll.setText('Add all')
            # QApplication.processEvents()
        except:
            self.addAll.setStyleSheet('color: red')
            self.addAll.setText('Error')
            pass

    def align_min_pdb_chain(self):
        selected = self.alignpdbChainAdded.selectedItems()
        dic = dict(self.alignaddedPdbChainDic)
        for k, v in dic.items():
            if v[0] in [i.text() for i in selected]:
                del self.alignaddedPdbChainDic[k]
        self.alignpdbChainAdded.clear()
        self.alignpdbChainAdded.addItems([v[0] for v in self.alignaddedPdbChainDic.values()])

    def check_local(self):
        # check local & similarity_flag
        if self.localGlobal.currentText() == 'Yes' and self.alignType.currentText() == 'Pairwise':
            self.similarityFlag.setCurrentIndex(0)
        elif self.localGlobal.currentText() == 'No' and self.alignType.currentText() == 'Pairwise':
            pass

        else:
            self.similarityFlag.setCurrentIndex(1)
        # check for local alignment & matrix offset
        if self.localGlobal.currentText() == 'Yes':
            self.matrixOffsetlabel.setEnabled(True)
            self.matrixOffset.setEnabled(True)
        else:
            self.matrixOffsetlabel.setEnabled(False)
            self.matrixOffset.setEnabled(False)

    def check_auto_overhang(self):
        if self.autoOverhang.currentText() == 'Yes':
            self.autoOverHangLimlabel.setEnabled(True)
            self.autoOverhangLim.setEnabled(True)
            self.overHangFaclabel.setEnabled(True)
            self.overhangFactor.setEnabled(True)
        else:
            self.autoOverHangLimlabel.setEnabled(False)
            self.autoOverhangLim.setEnabled(False)
            self.overHangFaclabel.setEnabled(False)
            self.overhangFactor.setEnabled(False)

    def user_matrix(self):
        self.user_matrix_file = QFileDialog.getOpenFileName()[0]
        if self.user_matrix_file:
            self.userMatrixTypelabel.setEnabled(True)
            self.userMatrixType.setEnabled(True)
        else:
            self.userMatrixTypelabel.setEnabled(False)
            self.userMatrixType.setEnabled(False)

    def align_get_values_th(self):
        t = threading.Thread(target=self.align_get_values)
        t.start()
        self.alignBut.setEnabled(False)

    def align_get_values(self, *args):
        try:
            # local alignment
            if self.localGlobal.currentText() == 'Yes':
                self.getLocalAlignL.insert(0, True)
                self.getLocalAlignL = self.getLocalAlignL[:1]
            else:
                self.getLocalAlignL.insert(0, False)
                self.getLocalAlignL = self.getLocalAlignL[:1]
            # alignment type >> self.alignType.currentText().upper()
            # matrix
            for i in self.matrixNameL:
                if self.matrix.currentText() == i:
                    self.getMatrixL.insert(0, self.matrixFileL[self.matrixNameL.index(i)])
                    self.getMatrixL = self.getMatrixL[:1]
            # gap function
            if self.gapFunction.currentText() == 'Yes':
                self.getGapFunctionL.insert(0, True)
                self.getGapFunctionL = self.getGapFunctionL[:1]
            else:
                self.getGapFunctionL.insert(0, False)
                self.getGapFunctionL = self.getGapFunctionL[:1]
            # gap open 1d >> float(self.gapOpen1d.text())
            # gap extension 1d >> float(self.gapExten1d.text())
            # gap open 3d >> float(self.gapOpen3d.text())
            # gap extension 3d >> float(self.gapExten3d.text())
            # gap residue score >> float(self.gapResScore.text())
            # gap gap score >> float(self.gapGapScore.text())
            # overhang >> int(self.overhang.text())
            # auto overhang
            if self.autoOverhang.currentText() == 'Yes':
                self.getAutoOverhangL.insert(0, True)
                self.getAutoOverhangL = self.getAutoOverhangL[:1]
            else:
                self.getAutoOverhangL.insert(0, False)
                self.getAutoOverhangL = self.getAutoOverhangL[:1]
            # auto overhang limit >> int(self.autoOverhangLim.text())
            # overhang factor >> float(self.overhangFactor.text())
            # similarity flag
            if self.similarityFlag.currentText() == 'Yes':
                self.getSimilarityFlagL.insert(0, True)
                self.getSimilarityFlagL = self.getSimilarityFlagL[:1]
            else:
                self.getSimilarityFlagL.insert(0, False)
                self.getSimilarityFlagL = self.getSimilarityFlagL[:1]
            # rms cutoff >> float(self.rmsCutOff.text())
            # dendrogramm
            if self.dendrogram.currentText() == 'No':
                self.getDendrogramL.insert(0, '')
                self.getDendrogramL = self.getDendrogramL[:1]
            else:
                fileName = os.path.join(self.path, 'dendrogramm_file.tree')
                self.getDendrogramL = self.getDendrogramL[:0]
                self.getDendrogramL.insert(0, fileName)
            # user matrix
            if self.user_matrix_file != '':
                self.getUsermatrixL.insert(0, self.user_matrix_file)
                self.getUsermatrixL = self.getUsermatrixL[:1]
            else:
                self.getUsermatrixL.insert(0, '')
                self.getUsermatrixL = self.getUsermatrixL[:1]
            # user matrix type >> self.userMatrixType.text().upper()
            # use fit
            if self.useFit.currentText() == 'No':
                self.getUseFitL.insert(0, False)
                self.getUseFitL = self.getUseFitL[:1]
            else:
                self.getUseFitL.insert(0, True)
                self.getUseFitL = self.getUseFitL[:1]
            # write fit
            if self.writeFit.currentText() == 'No':
                self.getWriteFitL.insert(0, False)
                self.getWriteFitL = self.getWriteFitL[:1]
            else:
                self.getWriteFitL.insert(0, True)
                self.getWriteFitL = self.getWriteFitL[:1]
        except Exception as er:
            print(er)
        self.run_alignment()

    def remove_ali_file(self):
        try:
            for file in self.alignPDB.selectedItems():
                os.unlink(os.path.join(self.path, file.text()))
            self.alignPDB.clear()
            self.alignPDB.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        except:
            pass

    def run_alignment(self):
        self.msgLabel2.setStyleSheet('color: black')
        self.msgLabel2.setText('Processing...')
        if self.seq.currentText() == 'Select':
            target = os.path.join(self.path, self.infoFile['Target'])
        if self.seq.currentText() != 'Select':
            target = os.path.join(self.path, self.seq.currentText())
        if len(self.ali_file) > 0:
            target = self.ali_file
        modeller_path = config.Config.get_modeller_path()
        sys.path.insert(0, modeller_path)
        try:
            from modeller import log, environ, alignment, model
            log.verbose()
            env = environ()
            env.io.atom_files_directory = './:../atom_files/'

            aln = alignment(env)
            for v in self.alignaddedPdbChainDic.values():
                mdl = model(env, file=v[1], model_segment=('FIRST:' + v[2][0], 'LAST:' + v[2][len(v[2]) - 1]))
                aln.append_model(mdl, atom_files=os.path.relpath(v[1]), align_codes=os.path.basename(os.path.splitext(v[1])[0]) + ':' + ''.join([i for i in v[2]]))
            if self.alignType.currentText() != 'Pairwise':
                for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                                    ((1., 1., 1., 1., 1., 0.), False, False)):
                    aln.salign(similarity_flag=self.getSimilarityFlagL[0], local_alignment=self.getLocalAlignL[0],
                               rms_cutoff=float(self.rmsCutOff.text()), normalize_pp_scores=False,
                               rr_file=self.getMatrixL[0], overhang=int(self.overhang.text()), auto_overhang=self.getAutoOverhangL[0],
                               overhang_auto_limit=int(self.autoOverhangLim.text()), overhang_factor=float(self.overhangFactor.text()),
                               gap_function=self.getGapFunctionL[0],
                               gap_penalties_1d=(float(self.gapOpen1d.text()), float(self.gapExten1d.text())),
                               gap_penalties_3d=(float(self.gapOpen3d.text()), float(self.gapExten3d.text())),
                               gap_gap_score=float(self.gapGapScore.text()), gap_residue_score=float(self.gapResScore.text()),
                               matrix_offset=float(self.matrixOffset.text()),
                               dendrogram_file=self.getDendrogramL[0], alignment_type=self.alignType.currentText().upper(),
                               input_weights_file=self.getUsermatrixL[0], weights_type=self.userMatrixType.currentText().upper(),
                               feature_weights=weights, improve_alignment=True, fit=self.getUseFitL[0], write_fit=write_fit,
                               write_whole_pdb=whole, output='ALIGNMENT QUALITY')
            else:
                for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                                    ((1., 1., 1., 1., 1., 0.), False, False)):
                    aln.salign(similarity_flag=self.getSimilarityFlagL[0], local_alignment=self.getLocalAlignL[0],
                               rms_cutoff=float(self.rmsCutOff.text()), normalize_pp_scores=False,
                               rr_file=self.getMatrixL[0], overhang=int(self.overhang.text()),
                               auto_overhang=self.getAutoOverhangL[0],
                               overhang_auto_limit=int(self.autoOverhangLim.text()),
                               overhang_factor=float(self.overhangFactor.text()),
                               gap_function=self.getGapFunctionL[0],
                               gap_penalties_1d=(float(self.gapOpen1d.text()), float(self.gapExten1d.text())),
                               gap_penalties_3d=(float(self.gapOpen3d.text()), float(self.gapExten3d.text())),
                               gap_gap_score=float(self.gapGapScore.text()), gap_residue_score=float(self.gapResScore.text()),
                               matrix_offset=float(self.matrixOffset.text()),
                               dendrogram_file=self.getDendrogramL[0], alignment_type=self.alignType.currentText().upper(),
                               input_weights_file=self.getUsermatrixL[0], weights_type=self.userMatrixType.currentText().upper(),
                               feature_weights=(1., 0., 0., 0., 0., 0.), improve_alignment=True, fit=self.getUseFitL[0], write_fit=write_fit,
                               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

            aln.write(file=os.path.join(self.path, 'MSA_templates.pap'), alignment_format='PAP')
            aln.write(file=os.path.join(self.path, 'MSA_templates.ali'), alignment_format='PIR')
            # The number of equivalent positions at different RMS_CUTOFF values can be
            # computed by changing the RMS value and keeping all feature weights = 0
            aln.salign(rms_cutoff=(1.0), normalize_pp_scores=False, rr_file=self.getMatrixL[0],
                       overhang=int(self.overhang.text()),
                       gap_penalties_1d=(float(self.gapOpen1d.text()), float(self.gapExten1d.text())),
                       gap_penalties_3d=(float(self.gapOpen3d.text()), float(self.gapExten3d.text())),
                       gap_gap_score=float(self.gapGapScore.text()), gap_residue_score=float(self.gapResScore.text()),
                       dendrogram_file=self.getDendrogramL[0], matrix_offset=float(self.matrixOffset.text()),
                       alignment_type=self.alignType.currentText().upper(), feature_weights=[0] * 6,
                       improve_alignment=False, fit=False, write_fit=self.getWriteFitL[0],
                       write_whole_pdb=False, output='QUALITY')
            # alignment of templates msa & target
            aln.clear()
            aln.append(file=os.path.join(self.path, 'MSA_templates.ali'), align_codes='all')
            aln_block = len(aln)
            aln.append(file=target, align_codes=os.path.splitext(os.path.basename(target))[0])
            aln.salign(max_gap_length=int(self.maxgapLen.text()),
                       gap_function=True,  # to use structure-dependent gap penalty
                       alignment_type='PAIRWISE', align_block=aln_block,
                       feature_weights=(1., 0., 0., 0., 0., 0.), overhang=int(self.overhang.text()),
                       gap_penalties_1d=(float(self.gapOpen1d.text()), float(self.gapExten1d.text())),
                       gap_penalties_2d=(float(self.gapPenalt2dHelix.text()), float(self.gapPenalt2dBeta.text()),
                                         float(self.gapPenalt2dAccess.text()), float(self.gapPenalt2dStright.text()),
                                         float(self.gapPenalt2dCACA.text()), float(self.gapPenalt2dDstMin.text()),
                                         float(self.gapPenalt2dDstPower.text()), float(self.gapPenalt2dStrucProfile.text()), 0.0),
                       similarity_flag=True)

            aln.write(file=os.path.join(self.path, self.fileName.text().strip() + '.aln'), alignment_format='PIR')
            aln.write(file=os.path.join(self.path, self.fileName.text().strip() + '.pap'), alignment_format='PAP')

            self.msgLabel2.setStyleSheet('color: green')
            self.msgLabel2.setText('Finished')
            self.alignBut.setEnabled(True)
        except ModuleNotFoundError:
            self.msgLabel2.setStyleSheet('color: red')
            self.msgLabel2.setText('MODELLER not found')
            self.alignBut.setEnabled(True)
        except Exception as er:
            self.msgLabel2.setStyleSheet('color: red')
            self.msgLabel2.setText('Error')
            self.alignBut.setEnabled(True)
            pass

    def view_align_th(self):
        t = threading.Thread(target=self.view_align)
        t.start()
        self.viewAlign.setEnabled(False)
        self.msgLabel2.setText('')

    def view_align(self):
        self.v = view_align.SequencesViewer()
        self.v.show()


# def main():
#     app = QApplication(sys.argv)
#     form = AlignTemp()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()
