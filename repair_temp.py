import os
import sys
import re
import json
import shutil
import threading
import collections
from Bio.PDB.PDBParser import PDBParser
from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog
import tpl_repair_temp
import config


class RepairTemp(QMainWindow, tpl_repair_temp.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    pdb_fileL = []
    aln_fileL = []
    YesNoL = ['Yes', 'No']
    patchTypeL = [['1MC', '+'], ['1MT', '+'], ['25P1', '+'], ['25P2', '+'], ['3PHO', '+'], ['3TER', '+'], ['5DP', '+'],
                ['5MC1', '+'], ['5MC2', '+'], ['5MET', '+'], ['5PHO', '+'], ['5TER', '+'], ['9EG', '+'],
                ['9MA', '+'], ['9MG', '+'], ['ACE', 1], ['ACP', 1], ['ASPP', '+'], ['CT1', 1], ['CT2', 1],
                ['CT3', 1], ['CTER', 1], ['DEO1', '+'], ['DEO2', '+'], ['DISU', 2], ['FHEM', '+'], ['GLUP', '+'],
                ['GLYP', 1], ['INO1', '+'], ['LIG1', 2], ['LIG2', 2], ['LIG3', 2], ['LINK', 2], ['NTER', 1],
                ['PHEM', 2], ['PLIG', 3], ['PLO2', 3], ['PROP', 1], ['PURA', '+'], ['PURG', '+'],
                ['PYRC', '+'], ['PYRU', '+'], ['TP1', '+'], ['TP1A', '+'], ['TP2', '+'], ['TP2A', '+']]
    chainL = []
    addedPdbChainResNumDic = collections.OrderedDict()
    oneResL = []
    twoResL = []
    threeResL = []
    resNumberL = []
    reNumberResL = []
    segmentResL = []
    patchTerminalL = []

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        # pdb files
        self.PDB.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        self.PDB.clicked.connect(self.get_pdb_chain)
        # browse pdb files
        self.browsePDB.released.connect(self.browse_pdb_th)
        # renumber res
        self.reNumberRes.addItems([i for i in self.YesNoL])
        self.reNumberRes.setCurrentIndex(1)
        # segment res
        self.segmentRes.addItems([i for i in self.YesNoL])
        self.segmentRes.setCurrentIndex(1)
        self.segmentRes.currentIndexChanged.connect(self.segment_res)
        # hetatm
        self.removeHETATM.addItems([i for i in self.YesNoL])
        self.removeHETATM.setCurrentIndex(1)
        # add patch
        self.addPatch.addItems([i for i in self.YesNoL])
        self.addPatch.setCurrentIndex(1)
        self.addPatch.currentIndexChanged.connect(self.use_patch)
        # patch terminal res
        self.patchTerminal.addItems([i for i in self.YesNoL])
        self.patchTerminal.setCurrentIndex(1)
        # click chain list
        self.chainPDB.clicked.connect(self.get_res_num)
        # type in search res
        self.searchRes.textEdited.connect(self.search_res)
        # patch types list
        self.patchType.addItems(i[0] for i in self.patchTypeL)
        # add patch
        self.addBut.released.connect(self.add_patch)
        # min patch
        self.minBut.released.connect(self.min_patch)
        # repair
        self.repairBut.released.connect(self.get_values_th)
        # remove file
        self.rep_remove_file.released.connect(self.remove_rep_file)

    def enterEvent(self, QEvent):
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']
        self.path = path
        #
        pdb_files = []
        for i in os.listdir(self.path):
            if i.endswith('.pdb') and i not in pdb_files:
                pdb_files.append(i)
        while self.pdb_fileL != pdb_files:
            self.PDB.clear()
            self.PDB.addItems([os.path.basename(i) for i in os.listdir(path) if i.endswith('.pdb')])
            self.pdb_fileL = pdb_files

    def browse_pdb_th(self):
        threading.Thread(target=self.browse_pdb).start()

    def browse_pdb(self):
        try:
            user_pdb = QFileDialog.getOpenFileNames()[0]
            if user_pdb:
                for file in user_pdb:
                    shutil.copy(os.path.abspath(file), self.path)
                self.PDB.clear()
                self.PDB.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        except:
            pass

    def segment_res(self):
        if self.segmentRes.currentText() == 'Yes':
            self.segFromLabel.setEnabled(True)
            self.segFrom.setEnabled(True)
            self.chainFrom.setEnabled(True)
            self.segTo.setEnabled(True)
            self.chainTo.setEnabled(True)
            self.get_pdb_chain()
        else:
            self.segFromLabel.setEnabled(False)
            self.segFrom.setEnabled(False)
            self.chainFrom.setEnabled(False)
            self.segTo.setEnabled(False)
            self.chainTo.setEnabled(False)

    def use_patch(self):
        if self.addPatch.currentText() == 'Yes':
            self.frame_addPatch.setEnabled(True)
            self.get_pdb_chain()
            self.search_res()
        else:
            self.frame_addPatch.setEnabled(False)

    def get_pdb_chain(self):
        try:
            self.chainL.clear()
            self.resNumberL.clear()
            self.chainFrom.clear()
            self.chainTo.clear()
            self.chainPDB.clear()
            self.resNumber.clear()

            p = PDBParser()
            if self.segmentRes.currentText() == 'Yes':
                s = p.get_structure(self.PDB.currentItem().text(),
                                    os.path.join(self.path, self.PDB.currentItem().text()))
                for i in s.get_chains():
                    self.chainL.append(i.id)

                self.chainFrom.addItems(self.chainL)
                self.chainFrom.setCurrentIndex(0)
                self.chainTo.addItems(self.chainL)
                self.chainTo.setCurrentIndex(0)

            if self.addPatch.currentText() == 'Yes':
                s = p.get_structure(self.PDB.currentItem().text(),
                                    os.path.join(self.path, self.PDB.currentItem().text()))
                for i in s.get_chains():
                    self.chainL.append(i.id)

                self.chainPDB.addItems(self.chainL)
                self.chainPDB.setCurrentRow(0)
                # get residue numbers
                for i in s[0][self.chainPDB.currentItem().text()]:
                    if i.id[0] == ' ':
                        self.resNumberL.append([str(i.id[1]) + '-' + i.get_resname()])
                self.resNumber.addItems(i[0] for i in self.resNumberL)
            self.search_res()
        except Exception as er:
            #print(er)
            pass

    def get_res_num(self):
        try:
            self.resNumberL.clear()
            self.resNumber.clear()

            p = PDBParser()
            s = p.get_structure(self.PDB.currentItem().text(), os.path.join(self.path, self.PDB.currentItem().text()))
            for i in s[0][self.chainPDB.currentItem().text()]:
                if i.id[0] == ' ':
                    self.resNumberL.append([str(i.id[1]) + '-' + i.get_resname()])
            self.resNumber.addItems(i[0] for i in self.resNumberL)
            self.search_res()
        except Exception as er:
            print(er)

    def search_res(self):
        searchL = []
        regSearch = re.compile(str(self.searchRes.text()).upper().strip())
        for i in self.resNumberL:
            if str(self.searchRes.text()).upper().strip() in regSearch.findall(str(i[0])):
                searchL.append(i)
        self.resNumber.clear()
        self.resNumber.addItems(i[0] for i in searchL)

    def add_patch(self):
        # add residues num & chain to list
        try:
            L = [self.chainPDB.currentItem().text(), self.str_to_num(self.resNumber.currentItem().text())]
            for i in self.patchTypeL:
                # patch for one ore more than 3 residues
                if i[0] == self.patchType.currentItem().text():
                    if i[1] == '+' or i[1] == 1:
                        self.oneResL.append(L)
                        k = [self.PDB.currentItem().text()] + [i for i in self.oneResL] + [self.patchType.currentItem().text()]
                        self.addedPdbChainResNumDic[os.path.splitext(k[0])[0] + ', ' + str(k[1][1]) + ': ' + k[1][0] + ', (' + k[2] + ')'] = [
                            k[0], k[1], k[2]]
                        self.pdbChainAdded.clear()
                        self.pdbChainAdded.addItems([k for k in self.addedPdbChainResNumDic.keys()])
                        # clear list for next key
                        self.oneResL.clear()

                    # patch for 2 residues
                    if i[1] == 2:
                        self.twoResL.append(L)
                        k = [self.PDB.currentItem().text()] + [i for i in self.twoResL[:2]] + [self.patchType.currentItem().text()]
                        try:
                            self.addedPdbChainResNumDic[
                                os.path.splitext(k[0])[0] + ', ' + str(k[1][1]) + ': ' + k[1][0] + ', ' + str(k[2][1]) + ': ' + k[2][0] + ' (' +
                                k[3] + ')'] = [k[0], k[1], k[2], k[3]]
                            self.pdbChainAdded.clear()
                            self.pdbChainAdded.addItems([k for k in self.addedPdbChainResNumDic.keys()])
                        except:
                            self.pdbChainAdded.clear()
                            self.pdbChainAdded.addItems([k for k in self.addedPdbChainResNumDic.keys()] + [
                                k[0] + ', ' + k[1][0] + ':' + str(k[1][1]) + ' (' + k[2] + ')'])
                            pass
                        # clear list for next key
                        if len(self.twoResL) == 2:
                            self.twoResL.clear()

                    # patch for 3 residues threeResL
                    if i[1] == 3:
                        self.threeResL.append(L)
                        k = [self.PDB.currentItem().text()] + [i for i in self.threeResL[:3]] + [self.patchType.currentItem().text()]
                        try:
                            self.addedPdbChainResNumDic[
                                os.path.splitext(k[0])[0] + ', ' + str(k[1][1]) + ': ' + k[1][0] + ', ' + str(k[2][1]) + ': ' + k[2][
                                    0] + ', ' + str(k[3][1]) + ': ' + k[3][0] + ' (' + k[4] + ')'] = [k[0], k[1], k[2],
                                                                                        k[3], k[4]]
                            self.pdbChainAdded.clear()
                            self.pdbChainAdded.addItems([k for k in self.addedPdbChainResNumDic.keys()])
                        except:
                            try:
                                self.pdbChainAdded.clear()
                                self.pdbChainAdded.addItems([k for k in self.addedPdbChainResNumDic.keys()] + [
                                    k[0] + ', ' + k[1][0] + ':' + str(k[1][1]) + ' (' + k[2] + ')'])
                            except:
                                self.pdbChainAdded.clear()
                                self.pdbChainAdded.addItems([k for k in self.addedPdbChainResNumDic.keys()] + [
                                    k[0] + ', ' + str(k[1][1]) + ': ' + k[1][0] + ', ' + str(k[2][1]) + ': ' + k[2][
                                        0] + ' (' + k[3] + ')'])
                                pass
                        # clear list for next key
                        if len(self.threeResL) == 3:
                            self.threeResL.clear()
        except Exception as er:
            #print(er)
            pass

    def min_patch(self):
        delete = []
        selected = self.pdbChainAdded.selectedItems()
        for i in selected:
            id = i.text()
            try:
                delete.append(id)
                del self.addedPdbChainResNumDic[id]
            except:
                pass
        self.pdbChainAdded.clear()
        self.pdbChainAdded.addItems([k for k in self.addedPdbChainResNumDic.keys()])

    def str_to_num(self, value):
        reg = re.compile(r'(\d+)')
        src = reg.search(value)
        return int(src.group(1))

    def remove_rep_file(self):
        try:
            for file in self.PDB.selectedItems():
                os.unlink(os.path.join(self.path, file.text()))
            self.PDB.clear()
            self.PDB.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        except:
            pass

    def remove_ali_file(self):
        try:
            self.PDB.clear()
            self.PDB.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        except:
            pass

    def get_values_th(self):
        t = threading.Thread(target=self.get_values)
        t.start()
        self.repairBut.setEnabled(False)
        self.msgLabel.setText('')

    def get_values(self):
        # get renumber residue
        if self.reNumberRes.currentText() == 'No':
            self.reNumberResL.insert(0, True)
            self.reNumberResL = self.reNumberResL[:1]
        elif self.reNumberRes.currentText() == 'Yes':
            self.reNumberResL.insert(0, False)
            self.reNumberResL = self.reNumberResL[:1]
        # get segment of residue
        if self.segmentRes.currentText() == 'Yes':
            if len(self.segFrom.text()) > 0 and len(self.segTo.text()) > 0:
                self.segmentResL.insert(0, (
                self.segFrom.text() + ':' + self.chainFrom.currentText(), self.segTo.text() + ':' + self.chainTo.currentText()))
                self.segmentResL = self.segmentResL[:1]
            else:
                self.segmentResL.insert(0, None)
                self.segmentResL = self.segmentResL[:1]
        else:
            self.segmentResL.insert(0, None)
            self.segmentResL = self.segmentResL[:1]
        # get patch terminal
        if self.patchTerminal.currentText() == 'No':
            self.patchTerminalL.insert(0, False)
            self.patchTerminalL = self.patchTerminalL[:1]
        elif self.patchTerminal.currentText() == 'yes':
            self.patchTerminalL.insert(0, True)
            self.patchTerminalL = self.patchTerminalL[:1]
        # get add patch >> get from self.addedPdbChainResNumDic
        # run repair
        self.repair_missing_atom()

    def repair_missing_atom(self):
        modeller_path = config.Config.get_modeller_path()
        sys.path.insert(0, modeller_path)
        try:
            from modeller import environ
            from modeller.scripts import complete_pdb
            env = environ()
            env.io.atom_files_directory = ['../atom_files']
            env.libs.topology.read(file='$(LIB)/top_heav.lib')
            env.libs.parameters.read(file='$(LIB)/par.lib')
            # remove HETATM
            if self.removeHETATM.currentText() == 'Yes':
                env.io.hetatm = False
            else:
                env.io.hetatm = True
            if len(self.addedPdbChainResNumDic) == 0:
                for n, file in enumerate(self.PDB.selectedItems()):
                    self.msgLabel.setStyleSheet('color: black')
                    self.msgLabel.setText('Repairing %d/%d...' % (n + 1, len(self.PDB.selectedItems())))
                    self.mdl = complete_pdb(env, os.path.join(self.path, file.text()), special_patches=None,
                                            transfer_res_num=self.reNumberResL[0], model_segment=self.segmentResL[0],
                                            patch_default=self.patchTerminalL[0])
                    newFile = os.path.splitext(file.text())[0] + '_repaired.pdb'
                    self.mdl.write(file=os.path.join(self.path, newFile), model_format='PDB')
                    print('File \'' + file.text() + '\' repaired.')
            elif self.addPatch.currentText() == 'Yes' and len(self.addedPdbChainResNumDic) > 0:
                # special patch
                def addpatch(mdl):
                    for v in self.addedPdbChainResNumDic.values():
                        if len(v) == 3:
                            mdl.patch(residue_type=v[2], residues=(mdl.residues[str(v[1][1]) + ':' + v[1][0]]))
                        if len(v) == 4:
                            mdl.patch(residue_type=v[3], residues=(
                            mdl.residues[str(v[1][1]) + ':' + v[1][0]], mdl.residues[str(v[2][1]) + ':' + v[2][0]]))
                        if len(v) == 5:
                            mdl.patch(residue_type=v[4], residues=(
                            mdl.residues[str(v[1][1]) + ':' + v[1][0]], mdl.residues[str(v[2][1]) + ':' + v[2][0]],
                            mdl.residues[str(v[3][1]) + ':' + v[3][0]]))

                self.mdl = complete_pdb(env, os.path.join(self.path, self.PDB.currentItem().text()),
                                        special_patches=addpatch, transfer_res_num=self.reNumberResL[0],
                                        model_segment=self.segmentResL[0], patch_default=self.patchTerminalL[0])
                newFile = os.path.splitext(self.PDB.currentItem().text())[0] + '_repaired.pdb'
                self.mdl.write(file=os.path.join(self.path, newFile), model_format='PDB')
                print('File \'' + self.PDB.currentItem().text() + '\' repaired.')
            # miss - get pdb file
            self.PDB.clear()
            self.PDB.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])

            self.msgLabel.setStyleSheet('color: green')
            self.msgLabel.setText('Finished')
            QApplication.processEvents()
            self.repairBut.setEnabled(True)
        except ModuleNotFoundError:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('MODELLER not found')
            self.repairBut.setEnabled(True)
        except Exception as er:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Error')
            self.repairBut.setEnabled(True)
            print(er)


# def main():
#     app = QApplication(sys.argv)
#     form = RepairTemp()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()
