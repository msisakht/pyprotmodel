import os
import sys
import re
import json
import threading
from winreg import *
from Bio.PDB import PDBParser
from PyQt5 import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication
import tpl_def_restraint


class DefineRestraint(QMainWindow, tpl_def_restraint.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    addedResChainDic = {}
    restraintTypeL = ['Stereo', 'Bond', 'Angle', 'Improper', 'Dihedral', 'Sphere', 'Sphere14', 'LJ', 'LJ14', 'Coulomb',
                      'Coulomb14', 'Distance', 'User_Distance', 'Nonb_Pair_Spline', 'Phi-Psi_Binormal', 'Phi_Dihedral',
                      'Psi_Dihedral', 'Omega_Dihedral', 'Chi1_Dihedral', 'Chi2_Dihedral', 'Chi3_Dihedral',
                      'Chi4_Dihedral']
    YesNoL = ['Yes', 'No']
    pdfWeightMethodL = ['Local', 'Global']
    alphaTypesL = ['α helix', 'β strand']
    addedResAtmNumDic = {}
    resNumberL = []
    atmFeaturesL = ['Distance', 'Angle', 'Dihedral angle', 'Solvent access', 'Density', 'X coordinate', 'Y coordinate',
                    'Z coordinate']
    matFormL = ['Lower bound', 'Upper bound', 'Gaussian', 'Multiple Gaussian', 'Lennard-Jones', 'Coulomb', 'Cosine']
    physicalFeatL = ['Bond', 'Angle', 'Dihedral', 'Improper', 'Soft_Sphere', 'Lennard_Jones', 'Coulomb', 'H_Bond',
                     'CA_Distance', 'N_O_Distance', 'Phi_Dihedral', 'Psi_Dihedral', 'Omega_Dihedral', 'Chi1_Dihedral',
                     'Chi2_Dihedral', 'Chi3_Dihedral', 'Chi4_Dihedral', 'Disulfide_Distance', 'Disulfide_Angle',
                     'Disulfide_Dihedral', 'Lower_Distance', 'Upper_Distance', 'SD_MN_Distance', 'Chi5_Dihedral',
                     'Phi_Psi_Dihedral', 'SD_SD_Distance', 'XY_Distance', 'NMR_Distance', 'NMR_Distance2',
                     'Min_Distance',
                     'Nonbond_Spline', 'Accessibility', 'Density', 'Absposition', 'Dihedral_Diff', 'GBSA', 'EM_Density',
                     'SAXS', 'Symmetry']
    addedFeatMatDic = {}
    getIntersegmentL = []
    getConvert2SplinL = []

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.PDBFile.addItems([i for i in os.listdir(self.path) if i.endswith('pdb')])
        self.PDBFile.clicked.connect(self.res_from_to)
        # add residue range
        self.addBut.released.connect(self.add_residue)
        # min residue range
        self.minBut.released.connect(self.min_residue)
        # restraint type
        self.restraintType.addItems(self.restraintTypeL)
        self.restraintType.currentIndexChanged.connect(self.rest_type_alignment)
        # alignment
        self.alignment.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.aln')])
        # accessibility type
        self.accessibilType.setText('8')
        # convert to spline
        self.convet2Splin.addItems(self.YesNoL)
        # interval size
        self.intervalSize.setText('0.5')
        # min interval
        self.minInterval.setText('5')
        # spline range
        self.splineRange.setText('4.0')
        # pdf weight method
        self.pdfWeightMethod.addItems(self.pdfWeightMethodL)
        # pdf weight cut off
        self.pdfCutOff.setText('0.05')
        # restraint intersegment
        self.intersegment.addItems(self.YesNoL)
        # max distance
        self.maxDistance.setText('6.0')
        # file name
        self.fileName.setText('rsr_file')
        # alpha/beta
        self.alphaTypes.addItems(self.alphaTypesL)
        # residue-number click
        self.resNumlist.clicked.connect(self.load_atom)
        # search residue
        self.searchRes.textEdited.connect(self.search_res)
        # add secondary structure
        self.secAddBut.released.connect(self.add_sec)
        # min secondary structure
        self.secMinBut.released.connect(self.min_sec)
        # atom feature
        self.atmFeatures.addItems(self.atmFeaturesL)
        # display atom level features
        self.atm3AndLabel.hide()
        self.atm3NameTo.hide()
        self.atm3NumTo.hide()
        self.atm3Chain.hide()
        self.atm3Check.hide()
        #
        self.atm4AndLabel.hide()
        self.atm4NameTo.hide()
        self.atm4NumTo.hide()
        self.atm4Chain.hide()
        self.atm4Check.hide()
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 525, 31))
        self.atmFeatures.currentIndexChanged.connect(self.display_atm_feature)
        # atom level feature, mean std,...
        self.matForm.addItems(self.matFormL)
        self.physicalFeat.addItems(self.physicalFeatL)
        self.mean.setText('3.5')
        self.stdev.setText('0.1')
        self.weights.setText('0')
        self.len_jonsParamA.setText('0')
        self.len_jonsParamB.setText('0')
        self.coulombQ1.setText('0')
        self.coulombQ2.setText('0')
        self.cosineForce.setText('0')
        self.cosinePhase.setText('0')
        self.cosinePeriod.setText('0')
        # display atom level math form
        self.weightsLabel.hide()
        self.weights.hide()
        self.len_jonsParamALabel.hide()
        self.len_jonsParamA.hide()
        self.len_jonsParamBLabel.hide()
        self.len_jonsParamB.hide()
        self.coulombQ1Label.hide()
        self.coulombQ1.hide()
        self.coulombQ2Label.hide()
        self.coulombQ2.hide()
        self.cosineForceLabel.hide()
        self.cosineForce.hide()
        self.cosinePhaseLabel.hide()
        self.cosinePhase.hide()
        self.cosinePeriodLabel.hide()
        self.cosinePeriod.hide()
        self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 380, 450, 31))
        self.matForm.currentIndexChanged.connect(self.display_atm_math_form)
        # add atom level features
        self.addFeature.released.connect(self.add_feat_mat_form)
        # min atom level features
        self.minFeature.released.connect(self.min_feat_mat_form)
        # check for pdb and ali files
        self.pdb_fileL = [os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')]
        self.aln_fileL = [os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.aln')]
        # apply & remove buttons
        self.applyBut.released.connect(self.get_value_th)
        self.removeBut.released.connect(self.remove_pdb_file)

    def get_modeller_path(self):
        try:
            pymodel_key = OpenKey(HKEY_CURRENT_USER, r'SOFTWARE\PyModel', 0, KEY_READ)
            [pathVal, regtype] = (QueryValueEx(pymodel_key, 'MODELLER_PATH'))
            CloseKey(pymodel_key)
            return pathVal
        except:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('MODELLER not found')

    def enterEvent(self, QEvent):
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']
        self.path = path
        #
        pdb_files = []
        aln_files = []
        for i in os.listdir(self.path):
            if i.endswith('.pdb') and i not in pdb_files:
                pdb_files.append(i)
            if i.endswith('.aln') and i not in aln_files:
                aln_files.append(i)
        while self.pdb_fileL != pdb_files:
            self.PDBFile.clear()
            self.PDBFile.addItems([os.path.basename(i) for i in os.listdir(path) if i.endswith('.pdb')])
            self.pdb_fileL = pdb_files
        while self.aln_fileL != aln_files:
            self.alignment.clear()
            self.alignment.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.aln')])
            self.aln_fileL = aln_files

    def res_from_to(self):
        self.resFrom.clear()
        self.chainFrom.clear()
        self.resTo.clear()
        self.chainTo.clear()
        self.alphaChainFrom.clear()
        self.alphaChainTo.clear()
        self.resNumlist.clear()
        self.resNumberL.clear()
        try:
            p = PDBParser()
            s = p.get_structure(self.PDBFile.currentItem().text(), os.path.join(self.path, self.PDBFile.currentItem().text()))
            self.resFrom.setText([str(i.id[1]) for i in s.get_residues()][0])
            self.resTo.insert([str(i.id[1]) for i in s.get_residues()][len(list(s.get_residues())) - 1])
            self.chainFrom.addItems(i.id for i in s.get_chains())
            self.chainTo.addItems(i.id for i in s.get_chains())
            self.alphaChainFrom.addItems(i.id for i in s.get_chains())
            self.alphaChainTo.addItems(i.id for i in s.get_chains())
            self.chainFrom.setCurrentIndex(0)
            self.chainTo.setCurrentIndex(0)
            self.alphaChainFrom.setCurrentIndex(0)
            self.alphaChainTo.setCurrentIndex(0)
            # get residue list
            for i in s[0].get_residues():
                self.resNumberL.append([str(i.id[1]) + '-' + i.get_resname()])
            if len(self.resNumberL) > 0:
                self.resNumlist.addItems(i[0] for i in self.resNumberL)
                self.resNumlist.setCurrentRow(0)
                # load atoms & chains
                self.search_res()
                self.load_atom()
        except Exception as er:
            #print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            pass

    def add_residue(self):
        try:
            k = [self.PDBFile.currentItem().text(), self.chainFrom.currentText(), self.resFrom.text().strip(),
                 self.chainTo.currentText(), self.resTo.text().strip()]
            # add data to dic
            if len(k[2].strip()) > 0 and len(k[4].strip()) > 0:
                if len(self.addedResChainDic) == 0:
                    self.addedResChainDic[k[0] + ', ' + k[2] + ': ' + k[1] + ', ' + k[4] + ': ' + k[3]] = [k[0], k[1], k[2],
                                                                                                           k[3], k[4]]
                else:
                    if k[0] == [v[0] for v in self.addedResChainDic.values()][0]:
                        self.addedResChainDic[k[0] + ', ' + k[2] + ': ' + k[1] + ', ' + k[4] + ': ' + k[3]] = [k[0], k[1],
                                                                                                               k[2], k[3],
                                                                                                               k[4]]
                self.resNum.clear()
                self.resNum.addItems([k for k in self.addedResChainDic.keys()])
        except:
            #print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            pass

    def min_residue(self):
        try:
            selected = self.resNum.selectedItems()
            for i in selected:
                del self.addedResChainDic[i.text()]
            self.resNum.clear()
            self.resNum.addItems([k for k in self.addedResChainDic.keys()])
        except:
            pass

    def rest_type_alignment(self):
        if self.restraintType.currentText().upper() in ('CHI1_DIHEDRAL', 'CHI2_DIHEDRAL', 'CHI3_DIHEDRAL', 'CHI4_DIHEDRAL',
                                                'PHI_DIHEDRAL', 'PSI_DIHEDRAL', 'OMEGA_DIHEDRAL', 'PHI-PSI_BINORMAL'):
            self.alignmentLabel.setEnabled(True)
            self.alignment.setEnabled(True)
            self.alignment.clear()
            self.alignment.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.aln')])
        else:
            self.alignmentLabel.setEnabled(False)
            self.alignment.setEnabled(False)
            self.alignment.clear()
            self.alignment.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.aln')])
            self.alignment.setCurrentIndex(0)

    def load_atom(self):
        chainL = []
        atomsL = []
        try:
            # get residue num
            r = re.search(r'(\d+)\s*-\s*(\w+)', self.resNumlist.currentItem().text())
            num = r.group(1)
            # get chain & load atoms
            p = PDBParser()
            s = p.get_structure(self.PDBFile.currentItem().text(), os.path.join(self.path, self.PDBFile.currentItem().text()))
            for i in s[0].get_residues():
                if i.resname == r.group(2) and i.id[1] == int(r.group(1)):
                    for at in i.get_atoms():
                        chainL.append(at.parent.parent.id)
                        atomsL.append(at.id)
            if len(chainL) > 0:
                pass
            else:
                chainL = ['']
            if not self.betaShCheckFrom.isChecked():
                self.betaShNameFrom.clear()
                self.betaShNameFrom.addItems(atomsL)
                self.betaShNumFrom.clear()
                self.betaShNumFrom.setText(num)
                self.betaShChainFrom.clear()
                self.betaShChainFrom.setText(chainL[0])
            if not self.betaShCheckTo.isChecked():
                self.betaShNameTo.clear()
                self.betaShNameTo.addItems(atomsL)
                self.betaShNumTo.clear()
                self.betaShNumTo.setText(num)
                self.betaShChainTo.clear()
                self.betaShChainTo.setText(chainL[0])
            if not self.atm1Check.isChecked():
                self.atm1Namefrom.clear()
                self.atm1Namefrom.addItems(atomsL)
                self.atm1Numfrom.clear()
                self.atm1Numfrom.setText(num)
                self.atm1Chain.clear()
                self.atm1Chain.setText(chainL[0])
            if not self.atm2Check.isChecked():
                self.atm2NameTo.clear()
                self.atm2NameTo.addItems(atomsL)
                self.atm2NumTo.clear()
                self.atm2NumTo.setText(num)
                self.atm2Chain.clear()
                self.atm2Chain.setText(chainL[0])
            if not self.atm3Check.isChecked():
                self.atm3NameTo.clear()
                self.atm3NameTo.addItems(atomsL)
                self.atm3NumTo.clear()
                self.atm3NumTo.setText(num)
                self.atm3Chain.clear()
                self.atm3Chain.setText(chainL[0])
            if not self.atm4Check.isChecked():
                self.atm4NameTo.clear()
                self.atm4NameTo.addItems(atomsL)
                self.atm4NumTo.clear()
                self.atm4NumTo.setText(num)
                self.atm4Chain.clear()
                self.atm4Chain.setText(chainL[0])
        except Exception as er:
            self.betaShNameFrom.clear()
            self.betaShNameTo.clear()
            self.atm1Namefrom.clear()
            self.atm2NameTo.clear()
            self.atm3NameTo.clear()
            self.atm4NameTo.clear()
            self.betaShNameFrom.clear()
            self.betaShNameTo.clear()
            self.atm1Namefrom.clear()
            self.atm2NameTo.clear()
            self.atm3NameTo.clear()
            self.atm4NameTo.clear()
            self.betaShNumFrom.clear()
            self.betaShNumTo.clear()
            self.atm1Numfrom.clear()
            self.atm2NumTo.clear()
            self.atm3NumTo.clear()
            self.atm4NumTo.clear()
            self.betaShChainFrom.clear()
            self.betaShChainTo.clear()
            self.atm1Chain.clear()
            self.atm2Chain.clear()
            self.atm3Chain.clear()
            self.atm4Chain.clear()
            pass
    
    def search_res(self):
        searchL = []
        regSearch = re.compile(str(self.searchRes.text()).upper().strip())
        for i in self.resNumberL:
            if str(self.searchRes.text()).upper().strip() in regSearch.findall(str(i[0])):
                searchL.append(i)
        self.resNumlist.clear()
        self.resNumlist.addItems(i[0] for i in searchL)
        self.resNumlist.setCurrentRow(0)
        self.load_atom()

    def add_sec(self):
        try:
            k = [self.PDBFile.currentItem().text(), self.alphaChainFrom.currentText(), self.alphaResFrom.text().strip(),
                 self.alphaChainTo.currentText(), self.alphaResTo.text().strip(), self.alphaTypes.currentText(),
                 self.betaShNameFrom.currentText(), self.betaShNumFrom.text().strip(), self.betaShChainFrom.text().strip(),
                 self.betaShNameTo.currentText().strip(), self.betaShNumTo.text().strip(), self.betaShChainTo.text().strip(),
                 self.HBonds.text().strip()]
            # add data to dic
            if len(k[2].strip()) > 0 and len(k[4].strip()) > 0:
                self.addedResAtmNumDic['Res ' + k[2] + ': ' + k[1] + ' to ' + k[4] + ': ' + k[3] + ' (' + k[5] + ')'] = [
                    k[0], k[1], k[2], k[3], k[4], k[5]]
            self.addedResAtmNum.clear()
            self.addedResAtmNum.addItems([k for k in self.addedResAtmNumDic.keys()])

            if len(k[6].strip()) > 0 and len(k[7].strip()) > 0 and len(k[12].strip()) > 0:
                if int(k[12]) > 0:
                    bShtype = '(parallel β sheet)'
                else:
                    bShtype = '(anti parallel β sheet)'
                self.addedResAtmNumDic[
                    'Atm ' + k[6].upper() + ':' + k[7] + ':' + k[8].upper() + ' to ' + k[9].upper() + ':' + k[10] + ':' + k[
                        11].upper() + ' & ' + k[12] + ' H-bonds ' + bShtype] = [k[6], k[7], k[8], k[9], k[10], k[11],
                                                                                  k[12]]
            self.addedResAtmNum.clear()
            self.addedResAtmNum.addItems([k for k in self.addedResAtmNumDic.keys()])
        except:
            pass

    def min_sec(self):
        selected = self.addedResAtmNum.selectedItems()
        for i in selected:
            del self.addedResAtmNumDic[i.text()]
        self.addedResAtmNum.clear()
        self.addedResAtmNum.addItems([k for k in self.addedResAtmNumDic.keys()])

    def display_atm_feature(self):
        try:
            if self.atmFeatures.currentText() == 'Distance':
                self.labelAfterFeat.setText('in Å between:')
                # 1
                self.atmAndLabel.show()
                self.atm1Namefrom.show()
                self.atm1Numfrom.show()
                self.atm1Chain.show()
                self.atm1Check.show()
                # 2
                self.atm2NameTo.show()
                self.atm2NumTo.show()
                self.atm2Chain.show()
                self.atm2Check.show()
                # 3
                self.atm3AndLabel.hide()
                self.atm3NameTo.hide()
                self.atm3NumTo.hide()
                self.atm3Chain.hide()
                self.atm3Check.hide()
                # 4
                self.atm4AndLabel.hide()
                self.atm4NameTo.hide()
                self.atm4NumTo.hide()
                self.atm4Chain.hide()
                self.atm4Check.hide()
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 525, 31))
            elif self.atmFeatures.currentText() == 'Angle':
                self.labelAfterFeat.setText('in rad between:')
                # 1
                self.atmAndLabel.show()
                self.atm1Namefrom.show()
                self.atm1Numfrom.show()
                self.atm1Chain.show()
                self.atm1Check.show()
                # 2
                self.atm2NameTo.show()
                self.atm2NumTo.show()
                self.atm2Chain.show()
                self.atm2Check.show()
                # 3
                self.atm3AndLabel.show()
                self.atm3NameTo.show()
                self.atm3NumTo.show()
                self.atm3Chain.show()
                self.atm3Check.show()
                # 4
                self.atm4AndLabel.hide()
                self.atm4NameTo.hide()
                self.atm4NumTo.hide()
                self.atm4Chain.hide()
                self.atm4Check.hide()
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 689, 31))
            elif self.atmFeatures.currentText() == 'Dihedral angle':
                self.labelAfterFeat.setText('in rad between:')
                # 1
                self.atmAndLabel.show()
                self.atm1Namefrom.show()
                self.atm1Numfrom.show()
                self.atm1Chain.show()
                self.atm1Check.show()
                # 2
                self.atm2NameTo.show()
                self.atm2NumTo.show()
                self.atm2Chain.show()
                self.atm2Check.show()
                # 3
                self.atm3AndLabel.show()
                self.atm3NameTo.show()
                self.atm3NumTo.show()
                self.atm3Chain.show()
                self.atm3Check.show()
                # 4
                self.atm4AndLabel.show()
                self.atm4NameTo.show()
                self.atm4NumTo.show()
                self.atm4Chain.show()
                self.atm4Check.show()
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 916, 31))
            elif self.atmFeatures.currentText() == 'Solvent access':
                self.labelAfterFeat.setText('of area (Å2) exposed to solvent of:')
                # 1
                self.atmAndLabel.hide()
                self.atm1Namefrom.show()
                self.atm1Numfrom.show()
                self.atm1Chain.show()
                self.atm1Check.show()
                # 2
                self.atm2NameTo.hide()
                self.atm2NumTo.hide()
                self.atm2Chain.hide()
                self.atm2Check.hide()
                # 3
                self.atm3AndLabel.hide()
                self.atm3NameTo.hide()
                self.atm3NumTo.hide()
                self.atm3Chain.hide()
                self.atm3Check.hide()
                # 4
                self.atm4AndLabel.hide()
                self.atm4NameTo.hide()
                self.atm4NumTo.hide()
                self.atm4Chain.hide()
                self.atm4Check.hide()
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 420, 31))
            elif self.atmFeatures.currentText() == 'Density':
                self.labelAfterFeat.setText('to:')
                # 1
                self.atmAndLabel.hide()
                self.atm1Namefrom.show()
                self.atm1Numfrom.show()
                self.atm1Chain.show()
                self.atm1Check.show()
                # 2
                self.atm2NameTo.hide()
                self.atm2NumTo.hide()
                self.atm2Chain.hide()
                self.atm2Check.hide()
                # 3
                self.atm3AndLabel.hide()
                self.atm3NameTo.hide()
                self.atm3NumTo.hide()
                self.atm3Chain.hide()
                self.atm3Check.hide()
                # 4
                self.atm4AndLabel.hide()
                self.atm4NameTo.hide()
                self.atm4NumTo.hide()
                self.atm4Chain.hide()
                self.atm4Check.hide()
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 291, 31))
            elif self.atmFeatures.currentText() == 'X coordinate':
                self.labelAfterFeat.setText('in Å to:')
                # 1
                self.atmAndLabel.hide()
                self.atm1Namefrom.show()
                self.atm1Numfrom.show()
                self.atm1Chain.show()
                self.atm1Check.show()
                # 2
                self.atm2NameTo.hide()
                self.atm2NumTo.hide()
                self.atm2Chain.hide()
                self.atm2Check.hide()
                # 3
                self.atm3AndLabel.hide()
                self.atm3NameTo.hide()
                self.atm3NumTo.hide()
                self.atm3Chain.hide()
                self.atm3Check.hide()
                # 4
                self.atm4AndLabel.hide()
                self.atm4NameTo.hide()
                self.atm4NumTo.hide()
                self.atm4Chain.hide()
                self.atm4Check.hide()
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 310, 31))
            elif self.atmFeatures.currentText() == 'Y coordinate':
                self.labelAfterFeat.setText('in Å to:')
                # 1
                self.atmAndLabel.hide()
                self.atm1Namefrom.show()
                self.atm1Numfrom.show()
                self.atm1Chain.show()
                self.atm1Check.show()
                # 2
                self.atm2NameTo.hide()
                self.atm2NumTo.hide()
                self.atm2Chain.hide()
                self.atm2Check.hide()
                # 3
                self.atm3AndLabel.hide()
                self.atm3NameTo.hide()
                self.atm3NumTo.hide()
                self.atm3Chain.hide()
                self.atm3Check.hide()
                # 4
                self.atm4AndLabel.hide()
                self.atm4NameTo.hide()
                self.atm4NumTo.hide()
                self.atm4Chain.hide()
                self.atm4Check.hide()
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 310, 31))
            elif self.atmFeatures.currentText() == 'Z coordinate':
                self.labelAfterFeat.setText('in Å to:')
                # 1
                self.atmAndLabel.hide()
                self.atm1Namefrom.show()
                self.atm1Numfrom.show()
                self.atm1Chain.show()
                self.atm1Check.show()
                # 2
                self.atm2NameTo.hide()
                self.atm2NumTo.hide()
                self.atm2Chain.hide()
                self.atm2Check.hide()
                # 3
                self.atm3AndLabel.hide()
                self.atm3NameTo.hide()
                self.atm3NumTo.hide()
                self.atm3Chain.hide()
                self.atm3Check.hide()
                # 4
                self.atm4AndLabel.hide()
                self.atm4NameTo.hide()
                self.atm4NumTo.hide()
                self.atm4Chain.hide()
                self.atm4Check.hide()
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 340, 310, 31))
        except Exception as er:
            #print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            pass

    def display_atm_math_form(self):
        try:
            if self.matForm.currentText() == 'Lower bound':
                self.meanLabel.show()
                self.mean.show()
                self.stdevLabel.show()
                self.stdev.show()
                #
                self.weightsLabel.hide()
                self.weights.hide()
                #
                self.len_jonsParamALabel.hide()
                self.len_jonsParamA.hide()
                self.len_jonsParamBLabel.hide()
                self.len_jonsParamB.hide()
                #
                self.coulombQ1Label.hide()
                self.coulombQ1.hide()
                self.coulombQ2Label.hide()
                self.coulombQ2.hide()
                #
                self.cosineForceLabel.hide()
                self.cosineForce.hide()
                self.cosinePhaseLabel.hide()
                self.cosinePhase.hide()
                self.cosinePeriodLabel.hide()
                self.cosinePeriod.hide()
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 380, 450, 31))
            elif self.matForm.currentText() == 'Upper bound':
                self.meanLabel.show()
                self.mean.show()
                self.stdevLabel.show()
                self.stdev.show()
                #
                self.weightsLabel.hide()
                self.weights.hide()
                #
                self.len_jonsParamALabel.hide()
                self.len_jonsParamA.hide()
                self.len_jonsParamBLabel.hide()
                self.len_jonsParamB.hide()
                #
                self.coulombQ1Label.hide()
                self.coulombQ1.hide()
                self.coulombQ2Label.hide()
                self.coulombQ2.hide()
                #
                self.cosineForceLabel.hide()
                self.cosineForce.hide()
                self.cosinePhaseLabel.hide()
                self.cosinePhase.hide()
                self.cosinePeriodLabel.hide()
                self.cosinePeriod.hide()
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 380, 450, 31))
            elif self.matForm.currentText() == 'Gaussian':
                self.meanLabel.show()
                self.mean.show()
                self.stdevLabel.show()
                self.stdev.show()
                #
                self.weightsLabel.hide()
                self.weights.hide()
                #
                self.len_jonsParamALabel.hide()
                self.len_jonsParamA.hide()
                self.len_jonsParamBLabel.hide()
                self.len_jonsParamB.hide()
                #
                self.coulombQ1Label.hide()
                self.coulombQ1.hide()
                self.coulombQ2Label.hide()
                self.coulombQ2.hide()
                #
                self.cosineForceLabel.hide()
                self.cosineForce.hide()
                self.cosinePhaseLabel.hide()
                self.cosinePhase.hide()
                self.cosinePeriodLabel.hide()
                self.cosinePeriod.hide()
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 380, 450, 31))
            elif self.matForm.currentText() == 'Multiple Gaussian':
                self.meanLabel.show()
                self.mean.show()
                self.stdevLabel.show()
                self.stdev.show()
                #
                self.weightsLabel.show()
                self.weights.show()
                #
                self.len_jonsParamALabel.hide()
                self.len_jonsParamA.hide()
                self.len_jonsParamBLabel.hide()
                self.len_jonsParamB.hide()
                #
                self.coulombQ1Label.hide()
                self.coulombQ1.hide()
                self.coulombQ2Label.hide()
                self.coulombQ2.hide()
                #
                self.cosineForceLabel.hide()
                self.cosineForce.hide()
                self.cosinePhaseLabel.hide()
                self.cosinePhase.hide()
                self.cosinePeriodLabel.hide()
                self.cosinePeriod.hide()
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 380, 530, 31))
            elif self.matForm.currentText() == 'Lennard-Jones':
                self.meanLabel.hide()
                self.mean.hide()
                self.stdevLabel.hide()
                self.stdev.hide()
                #
                self.weightsLabel.hide()
                self.weights.hide()
                #
                self.len_jonsParamALabel.show()
                self.len_jonsParamA.show()
                self.len_jonsParamBLabel.show()
                self.len_jonsParamB.show()
                #
                self.coulombQ1Label.hide()
                self.coulombQ1.hide()
                self.coulombQ2Label.hide()
                self.coulombQ2.hide()
                #
                self.cosineForceLabel.hide()
                self.cosineForce.hide()
                self.cosinePhaseLabel.hide()
                self.cosinePhase.hide()
                self.cosinePeriodLabel.hide()
                self.cosinePeriod.hide()
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 380, 490, 31))
            elif self.matForm.currentText() == 'Coulomb':
                self.meanLabel.hide()
                self.mean.hide()
                self.stdevLabel.hide()
                self.stdev.hide()
                #
                self.weightsLabel.hide()
                self.weights.hide()
                #
                self.len_jonsParamALabel.hide()
                self.len_jonsParamA.hide()
                self.len_jonsParamBLabel.hide()
                self.len_jonsParamB.hide()
                #
                self.coulombQ1Label.show()
                self.coulombQ1.show()
                self.coulombQ2Label.show()
                self.coulombQ2.show()
                #
                self.cosineForceLabel.hide()
                self.cosineForce.hide()
                self.cosinePhaseLabel.hide()
                self.cosinePhase.hide()
                self.cosinePeriodLabel.hide()
                self.cosinePeriod.hide()
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 380, 475, 31))
            elif self.matForm.currentText() == 'Cosine':
                self.meanLabel.hide()
                self.mean.hide()
                self.stdevLabel.hide()
                self.stdev.hide()
                #
                self.weightsLabel.hide()
                self.weights.hide()
                #
                self.len_jonsParamALabel.hide()
                self.len_jonsParamA.hide()
                self.len_jonsParamBLabel.hide()
                self.len_jonsParamB.hide()
                #
                self.coulombQ1Label.hide()
                self.coulombQ1.hide()
                self.coulombQ2Label.hide()
                self.coulombQ2.hide()
                #
                self.cosineForceLabel.show()
                self.cosineForce.show()
                self.cosinePhaseLabel.show()
                self.cosinePhase.show()
                self.cosinePeriodLabel.show()
                self.cosinePeriod.show()
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 380, 520, 31))
        except Exception as er:
            pass

    def add_feat_mat_form(self):
        try:
            k = [self.atmFeatures.currentText(), self.atm1Namefrom.currentText(), self.atm1Numfrom.text().strip(),
                    self.atm2NameTo.currentText(), self.atm2NumTo.text().strip(), self.atm3NameTo.currentText(),
                    self.atm3NumTo.text().strip(), self.atm4NameTo.currentText(), self.atm4NumTo.text().strip(),
                    self.matForm.currentText(), self.physicalFeat.currentText(), self.mean.text().strip(),
                    self.stdev.text().strip(), self.weights.text().strip(), self.len_jonsParamA.text().strip(),
                    self.len_jonsParamB.text().strip(), self.coulombQ1.text().strip(), self.coulombQ2.text().strip(),
                    self.cosinePhase.text().strip(), self.cosineForce.text().strip(), self.cosinePeriod.text().strip(),
                    self.atm1Chain.text().strip(), self.atm2Chain.text().strip(), self.atm3Chain.text().strip(),
                    self.atm4Chain.text().strip()]

            # add data to dic
            # ___ distance ___
            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[11], k[12], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[11], k[12], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[11], k[12], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + ' for each Gaussian. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' +
                        k[9] + ', ' + k[10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[11], k[12], k[13], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' +
                        k[9] + ', ' + k[10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[14], k[15], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + k[
                            9] + ', ' + k[10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[16], k[17], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For atoms: ' + k[1].upper() + ':' +
                        k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + ' (' + k[9] + ', ' + k[10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[18], k[19], k[20], k[21], k[22]]

            # ___ angle ___

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[11], k[12], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[11], k[12], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[11], k[12], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[11], k[12], k[13], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[14], k[15], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[16], k[17], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[18], k[19], k[20], k[21], k[22], k[23]]

            # ___ dihedral angle ___

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[11], k[12], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[11], k[12], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[11], k[12], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[11], k[12], k[13], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[14], k[15], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[16], k[17], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[18], k[19], k[20], k[21], k[22], k[23], k[24]]

            # ___ solvent access ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]
            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # ___ density ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]
            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # ___ X coordinate ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # ___ Y coordinate ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # ___ Z coordinate ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For atoms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # set to list box
            self.atmFeatMat.clear()
            self.atmFeatMat.addItems([k for k in self.addedFeatMatDic.keys()])
        except:
            pass

    def min_feat_mat_form(self):
        selected = self.atmFeatMat.selectedItems()
        for i in selected:
            del self.addedFeatMatDic[i.text()]
        self.atmFeatMat.clear()
        self.atmFeatMat.addItems([k for k in self.addedFeatMatDic.keys()])

    def get_value_th(self):
        t = threading.Thread(target=self.get_value)
        t.start()
        self.applyBut.setEnabled(False)
        self.msgLabel.setText('')

    def get_value(self):
        # restraint type >> self.restraintType.currentText()
        # alignment >> os.path.join(path, self.alignment.currentText())
        # intersegment
        if self.intersegment.currentText() == 'Yes':
           self.getIntersegmentL.insert(0, True)
           self.getIntersegmentL = self.getIntersegmentL[:1]
        elif self.intersegment.currentText() == 'No':
            self.getIntersegmentL.insert(0, False)
            self.getIntersegmentL = self.getIntersegmentL[:1]
        # atom range >> self.resRangeFrom.text() and self.resRangeTo.text()
        # accessibillity type >> int(self.accessibilType.currentText())
        # splin_on_site (convert 2 spline)
        if self.convet2Splin.currentText() == 'Yes':
           self.getConvert2SplinL.insert(0, True)
           self.getConvert2SplinL = self.getConvert2SplinL[:1]
        elif self.convet2Splin.currentText() == 'No':
            self.getConvert2SplinL.insert(0, False)
            self.getConvert2SplinL = self.getConvert2SplinL[:1]
        # interval size >> float(self.intervalSize.text())
        # min interval >> int(self.minInterval.text())
        # spline range >> float(self.splineRange.text())
        # pdf weight method >> self.pdfWeightMethod.currentText()
        # pdf cut off >> float(self.pdfCutOff.text())
        # max distance >> float(self.maxDistance.text())
        self.run_restraint()

    def remove_pdb_file(self):
        try:
            for file in self.PDBFile.selectedItems():
                os.unlink(os.path.join(self.path, file.text()))
            self.PDBFile.clear()
            self.PDBFile.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
            self.res_from_to()
            self.load_atom()
        except:
            pass

    def run_restraint(self):
        self.applyBut.setText('Processing...')
        modeller_path = self.get_modeller_path()
        sys.path.insert(0, modeller_path)
        try:
            from modeller import environ, alignment, log, selection, secondary_structure, physical, forms, features
            from modeller.scripts import complete_pdb
            res_ranges1 = [v for v in self.addedResChainDic.values()]
            res_ranges2 = [v for v in self.addedResAtmNumDic.values() if len(v) == 6]
            atm_ranges = [v for v in self.addedResAtmNumDic.values() if len(v) == 5]
            log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
            env = environ()
            env.io.atom_files_directory = '../atom_files'
            env.libs.topology.read(file='$(LIB)/top_heav.lib')
            env.libs.parameters.read(file='$(LIB)/par.lib')
            mdl = complete_pdb(env, os.path.join(self.path, self.PDBFile.currentItem().text()))
            aln = alignment(env)
            if self.restraintType.currentText().upper() in ('CHI1_DIHEDRAL', 'CHI2_DIHEDRAL', 'CHI3_DIHEDRAL','CHI4_DIHEDRAL', 'PHI_DIHEDRAL',
                                            'PSI_DIHEDRAL', 'OMEGA_DIHEDRAL', 'PHI-PSI_BINORMAL'):
                aln.append(file=os.path.join(self.path, self.alignment.currentText()))
            else:
                aln = None
            rsr = mdl.restraints
            s = selection()
            at = mdl.chains[0].atoms
            s.add(mdl.residue_range(num, chain) for num, chain in [(res_ranges1[n][2] + ':' + res_ranges1[n][1],
                            res_ranges1[n][4] + ':' + res_ranges1[n][3]) for n in range(len(res_ranges1))])
            rsr.make(atmsel=s, restraint_type=self.restraintType.currentText().upper(),
                     spline_on_site=self.getConvert2SplinL[0],
                     basis_pdf_weight=self.pdfWeightMethod.currentText().upper(),
                     basis_relative_weight=float(self.pdfCutOff.text()), intersegment=self.getIntersegmentL[0],
                     spline_dx=float(self.intervalSize.text()), spline_min_points=int(self.minInterval.text()),
                     spline_range=float(self.splineRange.text()), accessibility_type=int(self.accessibilType.text()),
                     distngh=float(self.maxDistance.text()), aln=aln)
            # add 2d structures
            for v in self.addedResAtmNumDic.values():
                if len(v) == 6 and v[5] == 'Alpha helix':
                    rsr.add(secondary_structure.alpha(mdl.residue_range(v[2] + ':' + v[1], v[4] + ':' + v[3])))
                if len(v) == 6 and v[5] == 'Beta strand':
                    rsr.add(secondary_structure.strand(mdl.residue_range(v[2] + ':' + v[1], v[4] + ':' + v[3])))
                if len(v) == 7:
                    rsr.add(secondary_structure.sheet(at[v[0].upper() + ':' + v[1] + ':' + v[2].upper()], at[v[3].upper() + ':' + v[4] + ':' + v[5].upper()], sheet_h_bonds=int(v[6])))
            # add atom restraints
            for v in self.addedFeatMatDic.values():
                # ___ get physical group ___
                if v[2] == 'Bond':
                    self.group = physical.bond
                if v[2] == 'Angle':
                    self.group = physical.angle
                if v[2] == 'Dihedral':
                    self.group = physical.dihedral
                if v[2] == 'Improper':
                    self.group = physical.improper
                if v[2] == 'Soft_Sphere':
                    self.group = physical.soft_sphere
                if v[2] == 'Lennard_Jones':
                    self.group = physical.lennard_jones
                if v[2] == 'Coulomb':
                    self.group = physical.coulomb
                if v[2] == 'H_Bond':
                    self.group = physical.h_bond
                if v[2] == 'CA_Distance':
                    self.group = physical.ca_distance
                if v[2] == 'N_O_Distance':
                    self.group = physical.n_o_distance
                if v[2] == 'Phi_Dihedral':
                    self.group = physical.phi_dihedral
                if v[2] == 'Psi_Dihedral':
                    self.group = physical.psi_dihedral
                if v[2] == 'Omega_Dihedral':
                    self.group = physical.omega_dihedral
                if v[2] == 'Chi1_Dihedral':
                    self.group = physical.chi1_dihedral
                if v[2] == 'Chi2_Dihedral':
                    self.group = physical.chi2_dihedral
                if v[2] == 'Chi3_Dihedral':
                    self.group = physical.chi3_dihedral
                if v[2] == 'Chi4_Dihedral':
                    self.group = physical.chi4_dihedral
                if v[2] == 'Disulfide_Distance':
                    self.group = physical.disulfide_distance
                if v[2] == 'Disulfide_Angle':
                    self.group = physical.disulfide_angle
                if v[2] == 'Disulfide_Dihedral':
                    self.group = physical.disulfide_dihedral
                if v[2] == 'Lower_Distance':
                    self.group = physical.lower_distance
                if v[2] == 'Upper_Distance':
                    self.group = physical.upper_distance
                if v[2] == 'SD_MN_Distance':
                    self.group = physical.sd_mn_distance
                if v[2] == 'Chi5_Dihedral':
                    self.group = physical.chi5_dihedral
                if v[2] == 'Phi_Psi_Dihedral':
                    self.group = physical.phi_psi_dihedral
                if v[2] == 'SD_SD_Distance':
                    self.group = physical.sd_sd_distance
                if v[2] == 'XY_Distance':
                    self.group = physical.xy_distance
                if v[2] == 'NMR_Distance':
                    self.group = physical.nmr_distance
                if v[2] == 'NMR_Distance2':
                    self.group = physical.nmr_distance2
                if v[2] == 'Min_Distance':
                    self.group = physical.min_distance
                if v[2] == 'Nonbond_Spline':
                    self.group = physical.nonbond_spline
                if v[2] == 'Accessibility':
                    self.group = physical.accessibility
                if v[2] == 'Density':
                    self.group = physical.density
                if v[2] == 'Absposition':
                    self.group = physical.absposition
                if v[2] == 'Dihedral_Diff':
                    self.group = physical.dihedral_diff
                if v[2] == 'GBSA':
                    self.group = physical.gbsa
                if v[2] == 'EM_Density':
                    self.group = physical.em_density
                if v[2] == 'SAXS':
                    self.group = physical.saxs
                if v[2] == 'Symmetry':
                    self.group = physical.symmetry

                # ___ distance ___
                if v[0] == 'Distance' and v[1] == 'Lower bound':
                    rsr.add(forms.lower_bound(group=self.group, feature=features.distance(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]), mean=float(v[7]),
                                              stdev=float(v[8])))

                if v[0] == 'Distance' and v[1] == 'Upper bound':
                    rsr.add(forms.upper_bound(group=self.group, feature=features.distance(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]), mean=float(v[7]),
                                              stdev=float(v[8])))

                if v[0] == 'Distance' and v[1] == 'Gaussian':
                    rsr.add(forms.gaussian(group=self.group, feature=features.distance(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]), mean=float(v[7]),
                                           stdev=float(v[8])))

                if v[0] == 'Distance' and v[1] == 'Multiple Gaussian':
                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.distance(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[10].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[11].upper()]), means=[float(v[7])],
                                                 stdevs=[float(v[8])], weights=[float(v[9])]))

                if v[0] == 'Distance' and v[1] == 'Lennard-Jones':
                    rsr.add(forms.lennard_jones(group=self.group, feature=features.distance(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]), A=float(v[7]), B=float(v[8])))

                if v[0] == 'Distance' and v[1] == 'Coulomb':
                    rsr.add(forms.coulomb(group=self.group, feature=features.distance(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]), q1=float(v[7]), q2=float(v[8])))

                if v[0] == 'Distance' and v[1] == 'Cosine':
                    rsr.add(forms.cosine(group=self.group, feature=features.distance(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[10].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]), phase=float(v[7]),
                                         force=float(v[8]), period=int(v[9])))

                # ___ angle ___

                if v[0] == 'Angle' and v[1] == 'Lower bound':
                    rsr.add(forms.lower_bound(group=self.group, feature=features.angle(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13].upper()]), mean=float(v[9]),
                                              stdev=float(v[10])))

                if v[0] == 'Angle' and v[1] == 'Upper bound':
                    rsr.add(forms.upper_bound(group=self.group, feature=features.angle(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13].upper()]), mean=float(v[9]),
                                              stdev=float(v[10])))

                if v[0] == 'Angle' and v[1] == 'Gaussian':
                    rsr.add(forms.gaussian(group=self.group, feature=features.angle(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13].upper()]), mean=float(v[9]),
                                           stdev=float(v[10])))

                if v[0] == 'Angle' and v[1] == 'Multiple Gaussian':
                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.angle(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[12].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[13].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[14].upper()]), means=[float(v[9])],
                                                 stdevs=[float(v[10])], weights=[float(v[11])]))

                if v[0] == 'Angle' and v[1] == 'Lennard-Jones':
                    rsr.add(forms.lennard_jones(group=self.group, feature=features.angle(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13].upper()]), A=float(v[9]), B=float(v[10])))

                if v[0] == 'Angle' and v[1] == 'Coulomb':
                    rsr.add(forms.coulomb(group=self.group,
                                          feature=features.angle(at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11]],
                                                                 at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12]],
                                                                 at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13]]),
                                          q1=float(v[9]), q2=float(v[10])))

                if v[0] == 'Angle' and v[1] == 'Cosine':
                    rsr.add(forms.cosine(group=self.group, feature=features.angle(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[12].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[13].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[14].upper()]), phase=float(v[9]),
                                         force=float(v[10]), period=int(v[11])))

                # ___ dihedral angle ___

                if v[0] == 'Dihedral angle' and v[1] == 'Lower bound':
                    rsr.add(forms.lower_bound(group=self.group, feature=features.dihedral(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[13].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[14].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[15].upper()],
                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]), mean=float(v[11]),
                                              stdev=float(v[12])))

                if v[0] == 'Dihedral angle' and v[1] == 'Upper bound':
                    rsr.add(forms.upper_bound(group=self.group, feature=features.dihedral(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[13].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[14].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[15].upper()],
                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]), mean=float(v[11]),
                                              stdev=float(v[12])))

                if v[0] == 'Dihedral angle' and v[1] == 'Gaussian':
                    rsr.add(forms.gaussian(group=self.group, feature=features.dihedral(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[13].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[14].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[15].upper()],
                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]), mean=float(v[11]),
                                           stdev=float(v[12])))

                if v[0] == 'Dihedral angle' and v[1] == 'Multiple Gaussian':
                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.dihedral(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[14].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[15].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[16].upper()],
                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[17].upper()]), means=[float(v[11])],
                                                 stdevs=[float(v[12])], weights=[float(v[13])]))

                if v[0] == 'Dihedral angle' and v[1] == 'Lennard-Jones':
                    rsr.add(forms.lennard_jones(group=self.group, feature=features.dihedral(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[13].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[14].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[15].upper()],
                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]), A=float(v[11]),
                                                B=float(v[12])))

                if v[0] == 'Dihedral angle' and v[1] == 'Coulomb':
                    rsr.add(forms.coulomb(group=self.group, feature=features.dihedral(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[13].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[14].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[15].upper()],
                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]), q1=float(v[11]),
                                          q2=float(v[12])))

                if v[0] == 'Dihedral angle' and v[1] == 'Cosine':
                    rsr.add(forms.cosine(group=self.group, feature=features.dihedral(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[14].upper()],
                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[15].upper()],
                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[16].upper()],
                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[17].upper()]), phase=float(v[11]),
                                         force=float(v[12]), period=int(v[13])))

                # ___ Solvent access ___

                if v[0] == 'Solvent access' and v[1] == 'Lower bound':
                    rsr.add(forms.lower_bound(group=self.group, feature=features.solvent_access(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Solvent access' and v[1] == 'Upper bound':
                    rsr.add(forms.upper_bound(group=self.group, feature=features.solvent_access(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Solvent access' and v[1] == 'Gaussian':
                    rsr.add(forms.gaussian(group=self.group, feature=features.solvent_access(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Solvent access' and v[1] == 'Multiple Gaussian':
                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.solvent_access(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), means=[float(v[5])],
                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                if v[0] == 'Solvent access' and v[1] == 'Lennard-Jones':
                    rsr.add(forms.lennard_jones(group=self.group, feature=features.solvent_access(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]), B=float(v[6])))

                if v[0] == 'Solvent access' and v[1] == 'Coulomb':
                    rsr.add(forms.coulomb(group=self.group, feature=features.solvent_access(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]), q2=float(v[6])))

                if v[0] == 'Solvent access' and v[1] == 'Cosine':
                    rsr.add(forms.cosine(group=self.group, feature=features.solvent_access(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), phase=float(v[5]), force=float(v[6]),
                                         period=int(v[7])))

                # ___ density ___

                if v[0] == 'Density' and v[1] == 'Lower bound':
                    rsr.add(forms.lower_bound(group=self.group, feature=features.density(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Density' and v[1] == 'Upper bound':
                    rsr.add(forms.upper_bound(group=self.group, feature=features.density(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Density' and v[1] == 'Gaussian':
                    rsr.add(forms.gaussian(group=self.group, feature=features.density(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Density' and v[1] == 'Multiple Gaussian':
                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.density(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), means=[float(v[5])],
                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                if v[0] == 'Density' and v[1] == 'Lennard-Jones':
                    rsr.add(forms.lennard_jones(group=self.group, feature=features.density(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]), B=float(v[6])))

                if v[0] == 'Density' and v[1] == 'Coulomb':
                    rsr.add(forms.coulomb(group=self.group, feature=features.density(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]), q2=float(v[6])))

                if v[0] == 'Density' and v[1] == 'Cosine':
                    rsr.add(forms.cosine(group=self.group,
                                         feature=features.density(at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                         phase=float(v[5]), force=float(v[6]), period=int(v[7])))

                # ___ x coordinate ___

                if v[0] == 'X coordinate' and v[1] == 'Lower bound':
                    rsr.add(forms.lower_bound(group=self.group, feature=features.x_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'X coordinate' and v[1] == 'Upper bound':
                    rsr.add(forms.upper_bound(group=self.group, feature=features.x_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'X coordinate' and v[1] == 'Gaussian':
                    rsr.add(forms.gaussian(group=self.group, feature=features.x_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'X coordinate' and v[1] == 'Multiple Gaussian':
                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.x_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), means=[float(v[5])],
                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                if v[0] == 'X coordinate' and v[1] == 'Lennard-Jones':
                    rsr.add(forms.lennard_jones(group=self.group, feature=features.x_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]), B=float(v[6])))

                if v[0] == 'X coordinate' and v[1] == 'Coulomb':
                    rsr.add(forms.coulomb(group=self.group, feature=features.x_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]), q2=float(v[6])))

                if v[0] == 'X coordinate' and v[1] == 'Cosine':
                    rsr.add(forms.cosine(group=self.group, feature=features.x_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), phase=float(v[5]), force=float(v[6]),
                                         period=int(v[7])))

                # ___ y coordinate ___

                if v[0] == 'Y coordinate' and v[1] == 'Lower bound':
                    rsr.add(forms.lower_bound(group=self.group, feature=features.y_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Y coordinate' and v[1] == 'Upper bound':
                    rsr.add(forms.upper_bound(group=self.group, feature=features.y_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Y coordinate' and v[1] == 'Gaussian':
                    rsr.add(forms.gaussian(group=self.group, feature=features.y_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Y coordinate' and v[1] == 'Multiple Gaussian':
                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.y_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), means=[float(v[5])],
                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                if v[0] == 'Y coordinate' and v[1] == 'Lennard-Jones':
                    rsr.add(forms.lennard_jones(group=self.group, feature=features.y_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]), B=float(v[6])))

                if v[0] == 'Y coordinate' and v[1] == 'Coulomb':
                    rsr.add(forms.coulomb(group=self.group, feature=features.y_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]), q2=float(v[6])))

                if v[0] == 'Y coordinate' and v[1] == 'Cosine':
                    rsr.add(forms.cosine(group=self.group, feature=features.y_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), phase=float(v[5]), force=float(v[6]),
                                         period=int(v[7])))

                # ___ z coordinate ___

                if v[0] == 'Z coordinate' and v[1] == 'Lower bound':
                    rsr.add(forms.lower_bound(group=self.group, feature=features.z_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Z coordinate' and v[1] == 'Upper bound':
                    rsr.add(forms.upper_bound(group=self.group, feature=features.z_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Z coordinate' and v[1] == 'Gaussian':
                    rsr.add(forms.gaussian(group=self.group, feature=features.z_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]), stdev=float(v[6])))

                if v[0] == 'Z coordinate' and v[1] == 'Multiple Gaussian':
                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.z_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), means=[float(v[5])],
                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                if v[0] == 'Z coordinate' and v[1] == 'Lennard-Jones':
                    rsr.add(forms.lennard_jones(group=self.group, feature=features.z_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]), B=float(v[6])))

                if v[0] == 'Z coordinate' and v[1] == 'Coulomb':
                    rsr.add(forms.coulomb(group=self.group, feature=features.z_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]), q2=float(v[6])))

                if v[0] == 'Z coordinate' and v[1] == 'Cosine':
                    rsr.add(forms.cosine(group=self.group, feature=features.z_coordinate(
                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]), phase=float(v[5]), force=float(v[6]),
                                         period=int(v[7])))

            if self.fileName.text().strip() == '':
                fileName = 'rsr_file'
            else:
                fileName = self.fileName.text().strip()
            rsr.write(os.path.join(self.path, fileName + '.rsr'))

            self.msgLabel.setStyleSheet('color: green')
            self.msgLabel.setText('Finished')
            self.applyBut.setText('Apply')
            self.applyBut.setEnabled(True)
        except ModuleNotFoundError:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('MODELLER not found')
            self.applyBut.setText('Apply')
            self.applyBut.setEnabled(True)
        except Exception as er:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Error')
            self.applyBut.setText('Apply')
            self.applyBut.setEnabled(True)

#
# def main():
#     app = QApplication(sys.argv)
#     form = DefineRestraint()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()