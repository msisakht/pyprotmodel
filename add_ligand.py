import os
import sys
import re
import json
import collections
from Bio.PDB import PDBParser
from Bio.SeqIO.PirIO import PirIterator
from PyQt5 import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication
import tpl_add_ligand


class AddLigand(QMainWindow, tpl_add_ligand.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    rsr_fileL = []
    ali_fileL = []
    non_temp_hetatmL = [['p', 'Pyroglutamate'], ['w', 'Water'], ['3', 'Calcium'], ['z', 'Zinc'], ['h', 'Heme'],
                        ['o', 'Oxygen'],
                        ['b', 'CO for heme'], ['c', 'Cystine'], ['H', 'Histidine'], ['&', 'ADP'], ['@', 'ATP'],
                        ['x', 'GDP'], ['y', 'GTP'], ['$', 'Magnesium'], ['d', 'Deoxyribose'], ['#', 'Dummy atom'],
                        ['0', 'Methylphosphate neutral'], ['1', 'Methylphosphate, anionic'],
                        ['2', 'Methylphosphate dianionic'], ['n', 'NAD'], ['4', 'NADH'], ['k', 'NIC'], ['5', 'NICH'],
                        ['f', 'PP'], ['r', 'Ribose'], ['i', 'SOD'], ['/', 'Chain break']]
    addedLigDic = collections.OrderedDict()
    non_Temp_addedLigDic = collections.OrderedDict()
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
    seqLeng = []
    residueL = []
    hetRes = []
    stdResL = []
    resLigCheck = ''

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.alignFile.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.aln')])
        self.alignFile.currentIndexChanged.connect(self.get_temps)
        # click on pdb
        self.templates.clicked.connect(self.get_std_lig)
        # click on std residue
        self.stdRes.clicked.connect(self.get_std_res_chain_atm)
        # click on lig
        self.tempLig.clicked.connect(self.get_lig_res_chain_atm)
        # non template ligands
        self.nTempLig.addItems(i[1] for i in self.non_temp_hetatmL)
        # add ligand
        self.addLigBut.released.connect(self.add_lig)
        # min ligand
        self.minLigBut.released.connect(self.min_lig)
        # hide restraint
        self.group_restraint.setEnabled(False)
        # add features
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
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 525, 31))
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
        self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 60, 450, 31))
        self.matForm.currentIndexChanged.connect(self.display_atm_math_form)
        # add atom level features
        self.addFeature_2.released.connect(self.add_feat_mat_form)
        # min atom level features
        self.minFeature_2.released.connect(self.min_feat_mat_form)
        # restraint file
        self.restraint.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.rsr')] + ['Define new'])
        self.restraint.currentIndexChanged.connect(self.def_restraint)
        # check for ali and rsr files
        self.ali_fileL = [i for i in os.listdir(self.path) if i.endswith('.aln')]
        self.rsr_fileL = [i for i in os.listdir(self.path) if i.endswith('.rsr')]
        # apply button
        self.applyBut.released.connect(self.add_lig_to_file)

    def enterEvent(self, QEvent):
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']
        self.path = path
        #
        ali_files = []
        rsr_files = []
        for i in os.listdir(self.path):
            if i.endswith('.aln') and i not in ali_files:
                ali_files.append(i)
            if i.endswith('.rsr') and i not in rsr_files:
                rsr_files.append(i)
        while self.ali_fileL != ali_files:
            self.alignFile.clear()
            self.alignFile.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.aln')])
            self.ali_fileL = ali_files
        while self.rsr_fileL != rsr_files:
            self.restraint.clear()
            self.restraint.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.rsr')] + ['Define new'])
            self.rsr_fileL = rsr_files

    def get_temps(self):
        self.templateDic = {}
        self.templates.clear()
        self.seqId.clear()
        templateL = []
        try:
            f = open(os.path.join(self.path, self.alignFile.currentText()))
            p = list(PirIterator(f))
            ids = [i.id for i in p]
            # set ids
            if self.alignFile.currentText() != 'Select':
                # add to dic
                for i in range(len(p)):
                    temp = p[i].description.split(':')[1]
                    self.templateDic[os.path.basename(temp)] = [os.path.basename(temp), p[i].id]
                    templateL.append(os.path.basename(temp))

                self.templates.addItems([i for i in templateL[:len(ids) - 1]])
                self.templates.setCurrentRow(0)
                # set sequence id
                self.seqId.setText(ids[len(ids) - 1:][0])
            self.get_std_lig()
        except Exception as er:
            #print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            self.templates.clear()
            pass

    def get_std_lig(self):
        self.resLigCheck = 'lig'
        self.hetResL = []
        self.stdRes.clear()
        self.tempLig.clear()
        resL = []
        ligL = []
        atomsL = []
        try:
            # get ligands, chains & atoms
            p = PDBParser()
            s = p.get_structure(self.templates.currentItem().text(), os.path.join(self.path, self.templates.currentItem().text()))
            for i in s[0].get_residues():
                if i.id[0] != ' ' and i.id[0] != 'W':
                    self.hetResL.append(i)
                    ligL.append(i.parent.id + ': ' + str(i.id[1]) + ' - ' + i.resname)
                    for at in i.get_atoms():
                        atomsL.append(at.id)
                if i.id[0] == ' ':
                    resL.append(i.parent.id + ': ' + str(i.id[1]) + ' - ' + i.resname)
            # check ligand
            self.stdRes.addItems(resL)
            self.tempLig.addItems(ligL)
            if len(ligL) > 0:
                self.tempLig.setCurrentRow(0)
            if len(ligL) == 0:
                self.tempLig.addItems('No Ligand')
                self.tempLig.setCurrentRow(0)
            if len(resL) > 0:
                self.stdRes.setCurrentRow(0)
            if len(resL) == 0:
                self.stdRes.addItems('No Residue')
                self.stdRes.setCurrentRow(0)
            # get residue num
            r = re.search(r'(\w*):\s*(\d+)\s+-\s+(\w+)', self.tempLig.currentItem().text())
            chain = r.group(1)
            num = r.group(2)
            liglenL = []
            for i in list(s[0].get_residues()):
                if i.id[0] != ' ' and i.id[0] != 'W':
                    liglenL.append(len(list(i.get_atoms())))
            atomsL = atomsL[:liglenL[0]]
            if self.group_restraint.isEnabled():
                if not self.atm1Check.isChecked():
                    self.atm1Namefrom.clear()
                    self.atm1Namefrom.addItems(atomsL)
                    self.atm1Numfrom.clear()
                    self.atm1Numfrom.setText(num)
                    self.atm1Chain.clear()
                    self.atm1Chain.setText(chain)
                if not self.atm2Check.isChecked():
                    self.atm2NameTo.clear()
                    self.atm2NameTo.addItems(atomsL)
                    self.atm2NumTo.clear()
                    self.atm2NumTo.setText(num)
                    self.atm2Chain.clear()
                    self.atm2Chain.setText(chain)
                if not self.atm3Check.isChecked():
                    self.atm3NameTo.clear()
                    self.atm3NameTo.addItems(atomsL)
                    self.atm3NumTo.clear()
                    self.atm3NumTo.setText(num)
                    self.atm3Chain.clear()
                    self.atm3Chain.setText(chain)
                if not self.atm4Check.isChecked():
                    self.atm4NameTo.clear()
                    self.atm4NameTo.addItems(atomsL)
                    self.atm4NumTo.clear()
                    self.atm4NumTo.setText(num)
                    self.atm4Chain.clear()
                    self.atm4Chain.setText(chain)
        except Exception as er:
            self.tempLig.clear()
            self.stdRes.clear()
            self.atm1Namefrom.clear()
            self.atm2NameTo.clear()
            self.atm3NameTo.clear()
            self.atm4NameTo.clear()
            #self.betaShNameFrom.clear()
            ##self.betaShNameTo.clear()
            self.atm1Namefrom.clear()
            self.atm2NameTo.clear()
            self.atm3NameTo.clear()
            self.atm4NameTo.clear()
            ##self.betaShNumFrom.clear()
            ##self.betaShNumTo.clear()
            self.atm1Numfrom.clear()
            self.atm2NumTo.clear()
            self.atm3NumTo.clear()
            self.atm4NumTo.clear()
            ##self.betaShChainFrom.clear()
            ##self.betaShChainTo.clear()
            self.atm1Chain.clear()
            self.atm2Chain.clear()
            self.atm3Chain.clear()
            self.atm4Chain.clear()
            self.stdRes.addItem('Error: %s' % (self.templates.currentItem().text()))
            pass

    def get_std_res_chain_atm(self):
        self.resLigCheck = 'stdRes'
        self.stdResL = []
        try:
            atomsL = []
            r = re.search(r'(\w*):\s*(\d+)\s+-\s+(\w+)', self.stdRes.currentItem().text())
            chain = r.group(1)
            num = r.group(2)
            res = r.group(3)

            # get ligands, chains & atoms
            p = PDBParser()
            s = p.get_structure(self.templates.currentItem().text(), os.path.join(self.path, self.templates.currentItem().text()))
            for i in s[0].get_residues():
                if i.parent.id == chain and str(i.id[1]) == num and i.resname == res:
                    self.stdResL.append(i)
                    for at in i.get_atoms():
                        atomsL.append(at.id)
            if self.group_restraint.isEnabled():
                if not self.atm1Check.isChecked():
                    self.atm1Namefrom.clear()
                    self.atm1Namefrom.addItems(atomsL)
                    self.atm1Numfrom.clear()
                    self.atm1Numfrom.setText(num)
                    self.atm1Chain.clear()
                    self.atm1Chain.setText(chain)
                if not self.atm2Check.isChecked():
                    self.atm2NameTo.clear()
                    self.atm2NameTo.addItems(atomsL)
                    self.atm2NumTo.clear()
                    self.atm2NumTo.setText(num)
                    self.atm2Chain.clear()
                    self.atm2Chain.setText(chain)
                if not self.atm3Check.isChecked():
                    self.atm3NameTo.clear()
                    self.atm3NameTo.addItems(atomsL)
                    self.atm3NumTo.clear()
                    self.atm3NumTo.setText(num)
                    self.atm3Chain.clear()
                    self.atm3Chain.setText(chain)
                if not self.atm4Check.isChecked():
                    self.atm4NameTo.clear()
                    self.atm4NameTo.addItems(atomsL)
                    self.atm4NumTo.clear()
                    self.atm4NumTo.setText(num)
                    self.atm4Chain.clear()
                    self.atm4Chain.setText(chain)
        except:
            self.tempLig.clear()
            if self.stdRes.currentItem().text().startswith('Error:'):
                pass
            else:
                self.stdRes.clear()
            self.atm1Namefrom.clear()
            self.atm2NameTo.clear()
            self.atm3NameTo.clear()
            self.atm4NameTo.clear()
            #self.betaShNameFrom.clear()
            #self.betaShNameTo.clear()
            self.atm1Namefrom.clear()
            self.atm2NameTo.clear()
            self.atm3NameTo.clear()
            self.atm4NameTo.clear()
            #self.betaShNumFrom.clear()
            #self.betaShNumTo.clear()
            self.atm1Numfrom.clear()
            self.atm2NumTo.clear()
            self.atm3NumTo.clear()
            self.atm4NumTo.clear()
            #self.betaShChainFrom.clear()
            #self.betaShChainTo.clear()
            self.atm1Chain.clear()
            self.atm2Chain.clear()
            self.atm3Chain.clear()
            self.atm4Chain.clear()
            pass

    def get_lig_res_chain_atm(self):
        self.resLigCheck = 'lig'
        atomsL = []
        try:
            p = PDBParser()
            s = p.get_structure(self.templates.currentItem().text(), os.path.join(self.path, self.templates.currentItem().text()))
            r = re.search(r'(\w):\s+(\d+)\s+-\s+\w+', self.tempLig.currentItem().text())
            chain = r.group(1)
            num = r.group(2)
            for i in s[0].get_residues():
                if str(i.id[1]) == r.group(2):
                    for at in i.get_atoms():
                        atomsL.append(at.id)
            if self.group_restraint.isEnabled():
                if not self.atm1Check.isChecked():
                    self.atm1Namefrom.clear()
                    self.atm1Namefrom.addItems(atomsL)
                    self.atm1Numfrom.clear()
                    self.atm1Numfrom.setText(num)
                    self.atm1Chain.clear()
                    self.atm1Chain.setText(chain)
                if not self.atm2Check.isChecked():
                    self.atm2NameTo.clear()
                    self.atm2NameTo.addItems(atomsL)
                    self.atm2NumTo.clear()
                    self.atm2NumTo.setText(num)
                    self.atm2Chain.clear()
                    self.atm2Chain.setText(chain)
                if not self.atm3Check.isChecked():
                    self.atm3NameTo.clear()
                    self.atm3NameTo.addItems(atomsL)
                    self.atm3NumTo.clear()
                    self.atm3NumTo.setText(num)
                    self.atm3Chain.clear()
                    self.atm3Chain.setText(chain)
                if not self.atm4Check.isChecked():
                    self.atm4NameTo.clear()
                    self.atm4NameTo.addItems(atomsL)
                    self.atm4NumTo.clear()
                    self.atm4NumTo.setText(num)
                    self.atm4Chain.clear()
                    self.atm4Chain.setText(chain)
        except:
            self.tempLig.clear()
            self.stdRes.clear()
            self.atm1Namefrom.clear()
            self.atm2NameTo.clear()
            self.atm3NameTo.clear()
            self.atm4NameTo.clear()
            #self.betaShNameFrom.clear()
            #self.betaShNameTo.clear()
            self.atm1Namefrom.clear()
            self.atm2NameTo.clear()
            self.atm3NameTo.clear()
            self.atm4NameTo.clear()
            #self.betaShNumFrom.clear()
            #self.betaShNumTo.clear()
            self.atm1Numfrom.clear()
            self.atm2NumTo.clear()
            self.atm3NumTo.clear()
            self.atm4NumTo.clear()
            #self.betaShChainFrom.clear()
            #self.betaShChainTo.clear()
            self.atm1Chain.clear()
            self.atm2Chain.clear()
            self.atm3Chain.clear()
            self.atm4Chain.clear()
            pass

    def add_lig(self):
        try:
            ligK = [self.templates.currentItem().text(), self.tempLig.currentItem().text()]
            try:
                if self.tempLig.currentItem().text() != '' and self.tempLig.currentItem().text() != 'No Ligand':
                    if len(self.addedLigDic) == 0:
                        self.addedLigDic[ligK[0] + ', ' + ligK[1]] = [ligK[0], os.path.splitext(ligK[0])[0], ligK[1]]
                    elif len(self.addedLigDic) > 0 and self.templates.currentItem().text() == [v[0] for v in self.addedLigDic.values()][0]:
                        self.addedLigDic[ligK[0] + ', ' + ligK[1]] = [ligK[0], os.path.splitext(ligK[0])[0], ligK[1]]
                    else:
                        pass
                # add non template ligands
                if self.nTempLig.currentItem().isSelected():
                    nLigK = [self.nTempLig.currentItem().text(), self.nTempLigCount.text().strip()]
                    if nLigK[1] == '':
                        nLigK[1] = '0'
                    self.non_Temp_addedLigDic['Non-template lig, ' + nLigK[0]] = [['Non-template lig, ' + nLigK[0] + ' (' + nLigK[1] + ')'], nLigK[0], nLigK[1]]
            except AttributeError:
                if len(self.addedLigDic) == 0:
                    self.addedLigDic[ligK[0] + ', ' + ligK[1]] = [ligK[0], os.path.splitext(ligK[0])[0], ligK[1]]
                elif len(self.addedLigDic) > 0 and self.templates.currentItem().text() == [v[0] for v in self.addedLigDic.values()][0]:
                    self.addedLigDic[ligK[0] + ', ' + ligK[1]] = [ligK[0], os.path.splitext(ligK[0])[0], ligK[1]]
                else:
                    pass
            self.ligandAdded.clear()
            self.ligandAdded.addItems([k for k in self.addedLigDic.keys()] + [v[0][0] for v in self.non_Temp_addedLigDic.values()])
        except Exception as er:
            #print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            pass

    def min_lig(self):
        selected = self.ligandAdded.selectedItems()
        for i in selected:
            if i.text() in self.addedLigDic.keys():
                del self.addedLigDic[i.text()]
            try:
                for k, v in self.non_Temp_addedLigDic.items():
                    if i.text() == v[0][0]:
                        del self.non_Temp_addedLigDic[k]
            except:
                continue
        self.ligandAdded.clear()
        self.ligandAdded.addItems([k for k in self.addedLigDic.keys()] + [v[0][0] for v in self.non_Temp_addedLigDic.values()])

    def def_restraint(self):
        if self.restraint.currentText() == 'Define new':
            self.group_restraint.setEnabled(True)
            self.get_lig_res_chain_atm()
        else:
            self.group_restraint.setEnabled(False)

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
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 525, 31))
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
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 689, 31))
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
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 916, 31))
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
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 420, 31))
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
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 291, 31))
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
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 310, 31))
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
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 310, 31))
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
                self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 310, 31))
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
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 60, 450, 31))
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
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 60, 450, 31))
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
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 60, 450, 31))
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
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 60, 530, 31))
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
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 60, 490, 31))
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
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 60, 475, 31))
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
                self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 60, 520, 31))
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
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[11], k[12], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[11], k[12], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[11], k[12], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + ' for each Gaussian. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' +
                        k[9] + ', ' + k[10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[11], k[12], k[13], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' +
                        k[9] + ', ' + k[10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[14], k[15], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + k[
                            9] + ', ' + k[10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[16], k[17], k[21], k[22]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '':
                if self.atmFeatures.currentText() == 'Distance' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For aroms: ' + k[1].upper() + ':' +
                        k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + str(k[22]).upper() + ' (' + ' (' + k[9] + ', ' + k[10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[18], k[19], k[20], k[21], k[22]]

            # ___ angle ___

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[11], k[12], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[11], k[12], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[11], k[12], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[11], k[12], k[13], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[14], k[15], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[16], k[17], k[21], k[22], k[23]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '':
                if self.atmFeatures.currentText() == 'Angle' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[18], k[19], k[20], k[21], k[22], k[23]]

            # ___ dihedral angle ___

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[11], k[12], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[11], k[12], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[11], k[12], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[11], k[12], k[13], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[14], k[15], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[16], k[17], k[21], k[22], k[23], k[24]]

            if k[1] != '' and k[2] != '' and k[3] != '' and k[4] != '' and k[5] != '' and k[6] != '' and k[7] != '' and k[8] != '':
                if self.atmFeatures.currentText() == 'Dihedral angle' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + str(k[21]).upper() + ', ' + k[3].upper() + ':' + k[4] + ':' + k[22].upper() + ', ' + k[5].upper() + ':' + k[6] + ':' + k[23].upper() + ', ' + k[7].upper() + ':' + k[8] + ':' + k[24].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[18], k[19], k[20], k[21], k[22], k[23], k[24]]

            # ___ solvent access ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]
            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Solvent access' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # ___ density ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]
            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Density' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # ___ X coordinate ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'X coordinate' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # ___ Y coordinate ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Y coordinate' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[18], k[19], k[20], k[21]]

            # ___ Z coordinate ___

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Lower bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be greater than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Upper bound':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be less than: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is harmonically restrained to be around: ' + k[11] + ' with standard deviation: ' +
                        k[12] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Multiple Gaussian':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a linear combination of Gaussians. ' + k[11] + ', ' + k[12] + ' and ' + k[
                            13] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[11], k[12], k[13], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Lennard-Jones':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of a Lennard-Jones potential, with control parameters A: ' + k[
                            14] + ' & B: ' + k[15] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[14], k[15], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Coulomb':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by means of an inverse square Coulomb potential created by charges: ' + k[
                            16] + ' & ' + k[17] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
                            10] + ').'] = [k[0], k[9], k[10], k[1], k[2], k[16], k[17], k[21]]

            if k[1] != '' and k[2] != '':
                if self.atmFeatures.currentText() == 'Z coordinate' and self.matForm.currentText() == 'Cosine':
                    self.addedFeatMatDic[
                        k[0] + ' is restrained by a CHARMM-style cosine function, with the phase shift: ' + k[
                            18] + ',force constant: ' + k[19] + ' & periodicity: ' + k[20] + '. For aroms: ' + k[1].upper() + ':' + k[2] + ':' + k[21].upper() + ' (' + k[9] + ', ' + k[
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

    def add_lig_to_file(self):
        try:
            if self.alignFile.currentText() != '':
                ligL = []
                gapL = []
                blkGapL = []
                blkL = []
                ligIndex = []
                length1 = []
                length2 = int()
                p = PDBParser(PERMISSIVE=1)
                s = p.get_structure(self.templates.currentItem().text(), os.path.join(self.path, self.templates.currentItem().text()))
                for i in s[0].get_residues():
                    if i.id[0] != ' ' and i.id[0] != 'W':
                        ligL.append(i.parent.id + ': ' + str(i.id[1]) + ' - ' + i.resname)
                for v in self.addedLigDic.values():
                    ligIndex.append(ligL.index(v[2]))
                for i in range(len(ligL)):
                    gapL.append('-')
                    blkL.append('.')
                for i in range(len(ligL)):
                    if i not in ligIndex:
                        blkGapL.append('-')
                    else:
                        blkGapL.append('.')
                # get data from alignment file
                f = open(os.path.join(self.path, self.alignFile.currentText()))
                p = list(PirIterator(f))
                pirCode = []
                for i in range(len(p)):
                    pirCode.append([v for v in p[i].annotations.values()])
                # pir data
                pirCode = [i[0] for i in pirCode]
                ids = [i.id for i in p]
                descriptions = [i.description for i in p]
                seqs = [i.seq for i in p]
                chainBreak = []
                if self.addChainIdent.isChecked():
                    chainBreak.append('/')
                else:
                    chainBreak.append('')
                # check for previously edited alignment
                for i in range(len(p)):
                    idx1 = p[i].description.split(':')[1]
                    reg = re.search(r'\((\+?\d+)\s*?,\s*?(\d+)\).*', p[i].description.split(':')[6])
                    if reg != None and self.templates.currentItem().text() == os.path.split(idx1)[1]:
                        length1 = [int(reg.group(1)), int(reg.group(2))]
                    elif reg != None and self.templates.currentItem().text() != os.path.split(idx1)[1] and len(length1) == 0:
                        length1 = [int(reg.group(2))]
                if len(length1) == 0:
                    length1.insert(0, len(seqs[0]))

                if len(length1) == 1:
                    length2 = length1[0]
                elif len(length1) > 1:
                    length2 = length1[1]

                # clear last added ligands
                with open(os.path.join(self.path, self.alignFile.currentText()), 'w') as f:
                    for i in range(len(p)):
                        idx1 = p[i].description.split(':')[1]
                        regClear = re.search(r'\((\+?\d+)\s*?,\s*?(\d+)\).*', p[i].description.split(':')[6])
                        if regClear != None and self.templates.currentItem().text() != os.path.split(idx1)[1]:
                            splt = p[i].description.split(':')
                            # get and update pdb length
                            idx4 = p[i].description.split(':')[4]
                            newIdx4 = regClear.group(1)
                            splt.remove(idx4)
                            splt.insert(4, newIdx4)
                            # new description
                            description2 = ':'.join([i for i in splt])

                            f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + description2 + '\n' + str(
                                seqs[i][:int(regClear.group(2))]) + '*' + '\n\n')
                        else:
                            try:
                                f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + str(
                                    seqs[i][:int(regClear.group(2))]) + '*' + '\n\n')
                            except:
                                pass

                # ___ add ligands to file ___

                # get non template ligands
                ntLigs = ''
                if len(self.non_Temp_addedLigDic) > 0:
                    ntLigsL = []
                    for v in self.non_Temp_addedLigDic.values():
                        for j in self.non_temp_hetatmL:
                            if v[1] == j[1]:
                                ntLigsL.append(j[0] * int(v[2]))
                    ntLigs = ''.join([i for i in ntLigsL])

                with open(os.path.join(self.path, self.alignFile.currentText()), 'w') as f:
                    try:
                        if len(self.addedLigDic) > 0 or len(self.non_Temp_addedLigDic) > 0:
                            for i in range(len(p)):
                                # get and add ligand for selected template
                                regaAdd = re.search(r'\((\+?\d+)\s*?,\s*?(\d+)\).*', p[i].description.split(':')[6])
                                if re.search(r'(.*):\w+|.*', ids[i]).group(1) in [v[1] for v in self.addedLigDic.values()]:
                                    splt = p[i].description.split(':')
                                    # get and update pdb length
                                    idx4 = p[i].description.split(':')[4]
                                    if len(length1) == 1:
                                        newIdx4 = str(int(idx4) + len(ligL))
                                        idx6 = splt[6]
                                        newIdx6 = '(' + idx4 + ', ' + str(length2) + ')' + splt[6]
                                        splt.remove(idx6)
                                        splt.insert(6, newIdx6)
                                    elif len(length1) > 1:
                                        newIdx4 = str(int(regaAdd.group(1)) + len(ligL))
                                    splt.remove(idx4)
                                    splt.insert(4, '+' + newIdx4)
                                    # new description
                                    description2 = ':'.join([i for i in splt])

                                    f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + description2 + '\n' + str(
                                        seqs[i][:length2]) + chainBreak[0] + ''.join(
                                        [i for i in blkL]) + '*' + '\n\n')
                                    print(ids[i], str(seqs[i][length2 - 10:length2]) + chainBreak[0] + ''.join(
                                        [i for i in blkL]) + '*')
                                # add siutable ligand to seqId seq
                                elif ids[i] == self.seqId.text().strip():
                                    f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + str(
                                        seqs[i][:length2]) + chainBreak[0] + ''.join(
                                        [i for i in blkGapL]) + ntLigs + '*' + '\n\n')
                                    print(ids[i], str(seqs[i][length2 - 10:length2]) + chainBreak[0] + ''.join(
                                        [i for i in blkGapL]) + ntLigs + '*\n')
                                # align other template with gaps
                                else:
                                    if regaAdd != None:
                                        splt = p[i].description.split(':')
                                        # get and update pdb length
                                        idx4 = p[i].description.split(':')[4]
                                        newIdx4 = regaAdd.group(1)
                                        splt.remove(idx4)
                                        splt.insert(4, newIdx4)
                                        # new description
                                        description2 = ':'.join([i for i in splt])

                                        f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + description2 + '\n' + str(
                                            seqs[i][:length2]) + chainBreak[0] + ''.join(
                                            [i for i in gapL]) + '*' + '\n\n')
                                        print(ids[i], str(seqs[i][length2 - 10:length2]) + chainBreak[0] + ''.join(
                                            [i for i in gapL]) + '*')

                                    else:
                                        f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + str(
                                            seqs[i][:length2]) + chainBreak[0] + ''.join(
                                            [i for i in gapL]) + '*' + '\n\n')
                                        print(ids[i], str(seqs[i][length2- 10:length2]) + chainBreak[0] + ''.join(
                                            [i for i in gapL]) + '*')
                            self.msgLabel.setStyleSheet('color: #000000')
                            self.msgLabel.setText('Ligand(s) added')
                        # clear all added ligands
                        elif len(self.addedLigDic) == 0 and len(self.non_Temp_addedLigDic) == 0:
                            for i in range(len(p)):
                                if len(length1) == 1:
                                    f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + str(
                                        seqs[i][:length2]) + '*' + '\n\n')
                                    print(ids[i], str(seqs[i][length2 - 10:length2]) + '*')

                                if len(length1) > 1:
                                    splt = p[i].description.split(':')
                                    # get and update pdb length
                                    idx1 = p[i].description.split(':')[1]
                                    idx4 = p[i].description.split(':')[4]
                                    newIdx4 = str(length1[0])
                                    splt.remove(idx4)
                                    splt.insert(4, '+' + newIdx4)
                                    # new description
                                    description2 = ':'.join([i for i in splt])
                                    if re.search(r'\((\+?\d+)\s*?,\s*?(\d+)\).*', p[i].description.split(':')[6]) != None and self.templates.currentItem().text() == os.path.split(idx1)[1]:
                                        f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + description2 + '\n' + str(
                                            seqs[i][:length2]) + '*' + '\n\n')
                                        print(ids[i], str(seqs[i][length2 - 10:length2]) + '*')
                                    else:
                                        f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + str(
                                            seqs[i][:length2]) + '*' + '\n\n')
                                        print(ids[i], str(seqs[i][length2 - 10:length2]) + '*')
                            print()
                            self.msgLabel.setStyleSheet('color: #000000')
                            self.msgLabel.setText('Alignment file restored')
                    except Exception as er:
                        self.msgLabel.setStyleSheet('color: red')
                        self.msgLabel.setText('Error')
                        print(er)

                # check for add ligands
                with open('info', 'r+') as f:
                    jFile = json.load(f)
                    if len(self.addedLigDic) > 0 or len(self.non_Temp_addedLigDic) > 0:
                        jFile['add_ligands'] = 'Yes'
                    elif len(self.addedLigDic) == 0 and len(self.non_Temp_addedLigDic) == 0:
                        jFile['add_ligands'] = 'No'
                        # add restraints
                    if self.restraint.currentText() != 'Select' and self.restraint.currentText() != 'Define new':
                        jFile['add_ligand_rsrFile'] = os.path.join(self.path, self.restraint.currentText())
                    elif self.restraint.currentText() == 'Select':
                        jFile['add_ligand_rsrFile'] = ''
                    if self.restraint.currentText() == 'Define new':
                        if len(self.addedFeatMatDic) > 0:
                            jFile['add_ligand_feature_mat'] = [v for v in self.addedFeatMatDic.values()]
                        else:
                            jFile['add_ligand_feature_mat'] = ''
                    elif self.restraint.currentText() == 'Select':
                        jFile['add_ligand_feature_mat'] = ''
                    f.seek(0)
                    f.write(json.dumps(jFile))
                    f.truncate()
        except Exception as er:
            #print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            pass


# def main():
#     app = QApplication(sys.argv)
#     form = AddLigand()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()