import os
import sys
import re
import json
import threading
from winreg import *
import pandas as pd
import shutil
from Bio.SeqIO.PirIO import PirIterator
import wget
from wget import bar_adaptive
from PyQt5 import QtGui
from PyQt5.QtWidgets import QMainWindow, QApplication, QHeaderView
import pd_model
import tpl_build_model


class BuildModel(QMainWindow, tpl_build_model.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    aln_fileL = []
    rsr_fileL = []
    vtfmL = ['Slow', 'Normal', 'Fast', 'Very fast', 'Fatest']
    assessL = ['GA341', 'DOPE', 'DOPE-HR', 'Normalized DOPE', 'SOAP-PP', 'SOAP-Loop', 'SOAP-Peptide', 'SOAP-Protein']
    YesNoL = ['Yes', 'No']
    alignFileL = []
    getVTFML = []
    chainsDic = {}
    addedResAtomDic = {}
    addedResAtomDicValues = []
    atomL = ['C', 'CA', 'CB', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2', 'CH2', 'CZ', 'CZ2',
             'CZ3', 'N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG',
             'OG1', 'OH', 'OXT', 'SD', 'SG']

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.alignFile.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.aln')])
        self.alignFile.currentIndexChanged.connect(self.get_temp)
        # deviation
        self.deviation.setText('4')
        # vtfm
        self.vtfm.addItems(self.vtfmL)
        self.vtfm.setCurrentIndex(1)
        # restraint file
        self.restraint.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.rsr')])
        # assess method
        self.assess.addItems(self.assessL)
        self.assess.clicked.connect(self.check_assess_method)
        # add hydrogen
        self.addHydrogen.addItems(self.YesNoL)
        self.addHydrogen.setCurrentIndex(1)
        # multi chain
        self.multiChain.addItems(self.YesNoL)
        self.multiChain.setCurrentIndex(1)
        self.multiChain.currentIndexChanged.connect(self.multi_chain)
        # align chain with model
        self.mChAlignWith.currentIndexChanged.connect(self.display_chain)
        # add/remove chains from aln file
        self.addchaintofileBut.released.connect(self.add_chain_to_file)
        self.minchainfromfileBut.released.connect(self.min_chain_to_file)
        # use symmetry
        self.frame_sym.setEnabled(False)
        self.useSymmetry.addItems(self.YesNoL)
        self.useSymmetry.setCurrentIndex(1)
        self.useSymmetry.currentIndexChanged.connect(self.use_sym)
        # symmetry atoms
        self.atmFrom.addItems(self.atomL)
        self.atmTo.addItems(self.atomL)
        # add/delete symmetry
        self.addsym_residue_atomBut.released.connect(self.add_sym)
        self.minsym_residue_atomBut.released.connect(self.min_sym)
        # number of models to build
        self.numModel.setText('5')
        # symmetry deviation
        self.symDev.setText('1')
        # check for ali and rsr files
        self.aln_fileL = [i for i in os.listdir(self.path) if i.endswith('.aln')]
        self.rsr_fileL = [i for i in os.listdir(self.path) if i.endswith('.rsr')]
        # save results as csv
        self.saveCsvBut.released.connect(self.save_result_csv)
        # hide result tableView
        self.group_model_result.hide()
        # download libraries
        self.downLibBut.hide()
        self.downLibBut.released.connect(self.download_libs)
        # build model & remove buttons
        self.buildBut.released.connect(self.build_model_th)

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

        while self.path != path:
            self.alignFile.clear()
            self.alignFile.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.aln')])
            #
            self.restraint.clear()
            self.restraint.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.rsr')])
            self.path = path
        #
        aln_files = []
        rsr_files = []
        for i in os.listdir(self.path):
            if i.endswith('.aln') and i not in aln_files:
                aln_files.append(i)
            if i.endswith('.rsr') and i not in rsr_files:
                rsr_files.append(i)
        while self.aln_fileL != aln_files:
            self.alignFile.clear()
            self.alignFile.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.aln')])
            self.aln_fileL = aln_files
        while self.rsr_fileL != rsr_files:
            self.restraint.clear()
            self.restraint.addItems(
                ['Select'] + [i for i in os.listdir(path) if i.endswith('.rsr')] + ['Define new'])
            self.rsr_fileL = rsr_files

    def get_temp(self):
        self.templates.clear()
        self.seqId.clear()
        self.chains.clear()
        self.chainFrom.clear()
        self.chainTo.clear()
        try:
            # get data from alignment file
            f = open(os.path.join(self.path, self.alignFile.currentText()))
            p = list(PirIterator(f))
            ids = [i.id for i in p]
            # set ids
            self.templates.addItems(ids[:len(ids) - 1])
            self.templates.setCurrentRow(0)
            self.seqId.setText(ids[len(ids) - 1])
            # check for multi chain templates & chains
            temps = [os.path.split(p[i].description.split(':')[1])[1] for i in range(len(p) - 1)]
            if self.multiChain.currentText() == 'Yes' and len(temps) > 0:
                self.mChAlignWith.clear()
                self.mChAlignWith.addItems(temps)
            else:
                self.mChAlignWith.clear()
            # set chain numbers
            try:
                for i in range(len(p) - 1):
                    chains = p[i].id
                    r = re.search(r'\w*:(\w*)', chains)
                    self.chainsDic.setdefault(os.path.split(p[i].description.split(':')[1])[1],
                                              ', '.join(i for i in r.group(1)))
            except:
                pass
            # load chains
            for k, v in self.chainsDic.items():
                if k == self.mChAlignWith.currentText():
                    c = re.findall(r'\w', v)
                    self.chains.setText(v)
                    self.chainFrom.addItems(c)
                    self.chainTo.addItems(c)
        except:
            pass

    def check_assess_method(self):
        soapL = [['soap_loop.hdf5', 'SOAP-Loop'], ['soap_peptide.hdf5', 'SOAP-Peptide'],
                 ['soap_protein_od.hdf5', 'SOAP-Protein']]
        no_soap = []
        for i in soapL:
            if os.path.exists(os.path.join(self.get_modeller_path(), i[0])):
                pass
            else:
                no_soap.append(i[1])
        for i in no_soap:
            if i == self.assess.currentItem().text():
                self.assess.currentItem().setForeground(QtGui.QColor('red'))
                self.assess.currentItem().setToolTip('Library not found.')
                self.downLibBut.show()

    def download_libs(self):
        soapL = [['soap_loop.hdf5', 'SOAP-Loop'], ['soap_peptide.hdf5', 'SOAP-Peptide'],
                 ['soap_protein_od.hdf5', 'SOAP-Protein']]
        for i in soapL:
            try:
                if i[1] in [i.text() for i in self.assess.selectedItems()]:
                    self.msgLabel.setStyleSheet('color: black')
                    self.msgLabel.setText('Downloading %s . . .' % (i[1]))
                    QApplication.processEvents()
                    wget.download('https://salilab.org/SOAP/' + i[0], self.get_modeller_path(), bar=bar_adaptive)
            except Exception as er:
                self.msgLabel.setStyleSheet('color: red')
                self.msgLabel.setText('Error')
                print(er)
                pass

    def multi_chain(self):
        if self.multiChain.currentText() == 'Yes':
            self.mChAlignWithLabel.setEnabled(True)
            self.mChAlignWith.setEnabled(True)
            self.chainsLabel.setEnabled(True)
            #self.chains.setEnabled(True)
            self.frame_addMinChain.setEnabled(True)
            self.get_temp()
        else:
            self.mChAlignWithLabel.setEnabled(False)
            self.mChAlignWith.setEnabled(False)
            self.chainsLabel.setEnabled(False)
            #self.chains.setEnabled(False)
            self.frame_addMinChain.setEnabled(False)

    def display_chain(self):
        self.chains.clear()
        self.chainFrom.clear()
        self.chainTo.clear()
        for k, v in self.chainsDic.items():
            if k == self.mChAlignWith.currentText():
                c = re.findall(r'\w', v)
                self.chains.setText(v)
                self.chainFrom.addItems(c)
                self.chainTo.addItems(c)

    def add_chain_to_file(self):
        try:
            f = open(os.path.join(self.path, self.alignFile.currentText()))
            p = list(PirIterator(f))
            mChSeq = [str(p[i].seq) for i in range(len(p)) if os.path.split(p[i].description.split(':')[1])[1] == self.mChAlignWith.currentText()][0]
            chainIdx = []
            # get index of '/'
            for i, j in enumerate(mChSeq):
                if j == '/':
                    chainIdx.append(i)
            # add '/' to seq
            l = [i for i in [str(p[i].seq) for i in range(len(p)) if p[i].id == self.seqId.text()][0]]
            if len(chainIdx) > 0:
                for i in chainIdx:
                    l.insert(i, '/')
            newSeq = ''.join([i for i in l])
            # add new seq to file
            pirCode = []
            for i in range(len(p)):
                pirCode.append([v for v in p[i].annotations.values()])
            # pir data
            pirCode = [i[0] for i in pirCode]
            ids = [i.id for i in p]
            descriptions = [i.description for i in p]
            seqs = [i.seq for i in p]
            with open(os.path.join(self.path, self.alignFile.currentText()), 'w') as f:
                for i in range(len(p)):
                    if ids[i] != self.seqId.text().strip():
                        f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + str(
                            seqs[i]) + '*' + '\n\n')
                        print(str(seqs[i]))
                    else:
                        f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + newSeq + '*' + '\n\n')

            self.msgLabel_chain.setStyleSheet('color: #000000')
            self.msgLabel_chain.setText('Chain(s) added')
        except:
            self.msgLabel_chain.setStyleSheet('color: red')
            self.msgLabel_chain.setText('Error')
            pass

    def min_chain_to_file(self):
        try:
            f = open(os.path.join(self.path, self.alignFile.currentText()))
            p = list(PirIterator(f))
            # remove '/' from seq
            pirCode = []
            for i in range(len(p)):
                pirCode.append([v for v in p[i].annotations.values()])
            # pir data
            pirCode = [i[0] for i in pirCode]
            ids = [i.id for i in p]
            descriptions = [i.description for i in p]
            seqs = [i.seq for i in p]
            newSeq = ''
            for i in range(len(p)):
                if p[i].id == self.seqId.text():
                    for char in str(p[i].seq):
                        if char != '/':
                            newSeq += char
            # write to file
            with open(os.path.join(self.path, self.alignFile.currentText()), 'w') as f:
                for i in range(len(p)):
                    if ids[i] != self.seqId.text().strip():
                        f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + str(
                            seqs[i]) + '*' + '\n\n')
                        print(str(seqs[i]))
                    else:
                        f.write('>' + pirCode[i] + ';' + ids[i] + '\n' + descriptions[i] + '\n' + newSeq + '*' + '\n\n')

            self.msgLabel_chain.setStyleSheet('color: #000000')
            self.msgLabel_chain.setText('Chain(s) removed')
        except:
            self.msgLabel_chain.setStyleSheet('color: red')
            self.msgLabel_chain.setText('Error')
            pass

    def use_sym(self):
        if self.useSymmetry.currentText() == 'Yes':
            self.frame_sym.setEnabled(True)
        else:
            self.frame_sym.setEnabled(False)

    def add_sym(self):
        try:
            k = [self.atmFrom.currentText(), self.chainFrom.currentText(), self.atmTo.currentText(), self.chainTo.currentText()]
            if len(k[0].strip()) > 0 and len(k[1].strip()) > 0 and len(k[2].strip()) > 0 and len(k[3].strip()) > 0:
                self.addedResAtomDic['(' + k[0].upper() + ': ' + k[1].upper() + ') | (' + k[2].upper() + ': ' + k[3].upper() + ')'] = [k[0], k[1], k[2], k[3], float(self.symDev.text().strip())]
            self.symAdded.clear()
            self.symAdded.addItems([k for k in self.addedResAtomDic.keys()])
            # add dic to json file
            with open('info', 'r+') as f:
                jFile = json.load(f)
                if len(self.addedResAtomDic) > 0:
                    jFile['symmetry'] = [v for v in self.addedResAtomDic.values()]
                else:
                    jFile['symmetry'] = []
                f.seek(0)
                f.write(json.dumps(jFile))
                f.truncate()
        except Exception as er:
            print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))

    def min_sym(self):
        try:
            selected = self.symAdded.selectedItems()
            for i in selected:
                del self.addedResAtomDic[i.text()]
            self.symAdded.clear()
            self.symAdded.addItems([k for k in self.addedResAtomDic.keys()])
            # update json file
            with open('info', 'r+') as f:
                jFile = json.load(f)
                if len(self.addedResAtomDic) > 0:
                    jFile['symmetry'] = [v for v in self.addedResAtomDic.values()]
                else:
                    jFile['symmetry'] = []
                f.seek(0)
                f.write(json.dumps(jFile))
                f.truncate()
        except:
            pass

    def save_result_csv(self):
        try:
            self.msglabel_csv.setStyleSheet('color: #000000')
            self.msglabel_csv.setText('Saved')
            QApplication.processEvents()
            pd.DataFrame.to_csv(self.model_df, os.path.join(self.path, 'models.csv'), index=False)
        except Exception as er:
            self.msglabel_csv.setStyleSheet('color: red')
            self.msglabel_csv.setText('Error')
            print(er)
            pass

    def build_model_th(self):
        t = threading.Thread(target=self.build_model)
        t.start()
        self.buildBut.setEnabled(False)
        self.msgLabel.setText('')

    def build_model(self):
        self.buildBut.setText('Processing...')
        modeller_path = self.get_modeller_path()
        sys.path.insert(0, modeller_path)
        try:
            from modeller import environ, soap_protein_od, soap_pp, soap_loop, soap_peptide
            from modeller.automodel import automodel, autosched, assess
            assessFuncL = (
            assess.GA341, assess.DOPE, assess.DOPEHR, assess.normalized_dope, soap_pp.Assessor(), soap_loop.Scorer(),
            soap_peptide.Scorer(), soap_protein_od.Scorer())

            # get align file
            self.alignFileL.insert(0, os.path.join(self.path, self.alignFile.currentText()))
            self.alignFileL = self.alignFileL[:1]
            # get temps >> self.templates.currentItem().text()
            # get seqId >> self.seqId.text()
            # get deviation >> int( self.deviation.text())
            # get vtfm (lib_schedule)
            lib_schedule = autosched.normal
            if self.vtfm.currentText() == 'Slow':
                lib_schedule = autosched.slow
            elif self.vtfm.currentText() == 'Normal':
                lib_schedule = autosched.normal
            elif self.vtfm.currentText() == 'Fast':
                lib_schedule = autosched.fast
            elif self.vtfm.currentText() == 'Very fast':
                lib_schedule = autosched.very_fast
            elif self.vtfm.currentText() == 'Fatest':
                lib_schedule = autosched.fastest
            # get assess method
            getAssessL = []
            selection = self.assess.selectedItems()
            for i in selection:
                getAssessL.append(i.text())
            idx = []
            for i in self.assessL:
                if i in getAssessL:
                    idx.append(self.assessL.index(i))
            assess_methods = ([i for i in assessFuncL if assessFuncL.index(i) in idx])
            # num of models >>  int(self.numModel.text())
            # check for move models to it's path
            fileL_before = [i for i in os.listdir(os.getcwd())]
            # ___ run automodel ____
            env = environ()
            if self.infoFile['add_ligands'] == 'Yes':
                env.io.hetatm = True
            else:
                env.io.hetatm = False
            if self.addHydrogen.currentText() == 'Yes':
                env.io.hydrogen = True
            else:
                env.io.hydrogen = False
            if self.restraint.currentText() == 'Select':
                restraint = None
            else:
                restraint = open(os.path.join(self.path, self.restraint.currentText()), 'r')

            class MyModel(automodel):

                def special_restraints(self, aln):
                    # get path
                    file = open('info')
                    infoFile = json.load(file)

                    sys.path.insert(0, modeller_path)
                    from modeller import physical, forms, features, symmetry, selection
                    rsr = self.restraints
                    at = self.atoms
                    try:
                        rsr.append(restraint)
                    except:
                        pass
                    if infoFile['add_ligand_rsrFile'] != '':
                        rsr.append(file=infoFile['add_ligand_rsrFile'])
                    elif infoFile['add_ligand_feature_mat'] != '':
                        for v in infoFile['add_ligand_feature_mat']:
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
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]),
                                                              mean=float(v[7]),
                                                              stdev=float(v[8])))

                                if v[0] == 'Distance' and v[1] == 'Upper bound':
                                    rsr.add(forms.upper_bound(group=self.group, feature=features.distance(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]),
                                                              mean=float(v[7]),
                                                              stdev=float(v[8])))

                                if v[0] == 'Distance' and v[1] == 'Gaussian':
                                    rsr.add(forms.gaussian(group=self.group, feature=features.distance(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]),
                                                           mean=float(v[7]),
                                                           stdev=float(v[8])))

                                if v[0] == 'Distance' and v[1] == 'Multiple Gaussian':
                                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.distance(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[10].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[11].upper()]),
                                                                 means=[float(v[7])],
                                                                 stdevs=[float(v[8])], weights=[float(v[9])]))

                                if v[0] == 'Distance' and v[1] == 'Lennard-Jones':
                                    rsr.add(forms.lennard_jones(group=self.group, feature=features.distance(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]), A=float(v[7]),
                                                                B=float(v[8])))

                                if v[0] == 'Distance' and v[1] == 'Coulomb':
                                    rsr.add(forms.coulomb(group=self.group, feature=features.distance(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[9].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]), q1=float(v[7]),
                                                          q2=float(v[8])))

                                if v[0] == 'Distance' and v[1] == 'Cosine':
                                    rsr.add(forms.cosine(group=self.group, feature=features.distance(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[10].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[10].upper()]),
                                                         phase=float(v[7]),
                                                         force=float(v[8]), period=int(v[9])))

                                # ___ angle ___

                                if v[0] == 'Angle' and v[1] == 'Lower bound':
                                    rsr.add(forms.lower_bound(group=self.group, feature=features.angle(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13].upper()]),
                                                              mean=float(v[9]),
                                                              stdev=float(v[10])))

                                if v[0] == 'Angle' and v[1] == 'Upper bound':
                                    rsr.add(forms.upper_bound(group=self.group, feature=features.angle(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13].upper()]),
                                                              mean=float(v[9]),
                                                              stdev=float(v[10])))

                                if v[0] == 'Angle' and v[1] == 'Gaussian':
                                    rsr.add(forms.gaussian(group=self.group, feature=features.angle(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13].upper()]),
                                                           mean=float(v[9]),
                                                           stdev=float(v[10])))

                                if v[0] == 'Angle' and v[1] == 'Multiple Gaussian':
                                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.angle(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[12].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[13].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[14].upper()]),
                                                                 means=[float(v[9])],
                                                                 stdevs=[float(v[10])], weights=[float(v[11])]))

                                if v[0] == 'Angle' and v[1] == 'Lennard-Jones':
                                    rsr.add(forms.lennard_jones(group=self.group, feature=features.angle(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13].upper()]), A=float(v[9]),
                                                                B=float(v[10])))

                                if v[0] == 'Angle' and v[1] == 'Coulomb':
                                    rsr.add(forms.coulomb(group=self.group, feature=features.angle(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[11]],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[12]],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[13]]), q1=float(v[9]),
                                                          q2=float(v[10])))

                                if v[0] == 'Angle' and v[1] == 'Cosine':
                                    rsr.add(forms.cosine(group=self.group, feature=features.angle(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[12].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[13].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[14].upper()]),
                                                         phase=float(v[9]),
                                                         force=float(v[10]), period=int(v[11])))

                                # ___ dihedral angle ___

                                if v[0] == 'Dihedral angle' and v[1] == 'Lower bound':
                                    rsr.add(forms.lower_bound(group=self.group, feature=features.dihedral(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[13].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[14].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[15].upper()],
                                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]),
                                                              mean=float(v[11]),
                                                              stdev=float(v[12])))

                                if v[0] == 'Dihedral angle' and v[1] == 'Upper bound':
                                    rsr.add(forms.upper_bound(group=self.group, feature=features.dihedral(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[13].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[14].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[15].upper()],
                                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]),
                                                              mean=float(v[11]),
                                                              stdev=float(v[12])))

                                if v[0] == 'Dihedral angle' and v[1] == 'Gaussian':
                                    rsr.add(forms.gaussian(group=self.group, feature=features.dihedral(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[13].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[14].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[15].upper()],
                                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]),
                                                           mean=float(v[11]),
                                                           stdev=float(v[12])))

                                if v[0] == 'Dihedral angle' and v[1] == 'Multiple Gaussian':
                                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.dihedral(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[14].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[15].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[16].upper()],
                                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[17].upper()]),
                                                                 means=[float(v[11])],
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
                                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[16].upper()]),
                                                          q1=float(v[11]),
                                                          q2=float(v[12])))

                                if v[0] == 'Dihedral angle' and v[1] == 'Cosine':
                                    rsr.add(forms.cosine(group=self.group, feature=features.dihedral(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[14].upper()],
                                        at[str(v[5]).upper() + ':' + str(v[6]) + ':' + v[15].upper()],
                                        at[str(v[7]).upper() + ':' + str(v[8]) + ':' + v[16].upper()],
                                        at[str(v[9]).upper() + ':' + str(v[10]) + ':' + v[17].upper()]),
                                                         phase=float(v[11]),
                                                         force=float(v[12]), period=int(v[13])))

                                # ___ Solvent access ___

                                if v[0] == 'Solvent access' and v[1] == 'Lower bound':
                                    rsr.add(forms.lower_bound(group=self.group, feature=features.solvent_access(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'Solvent access' and v[1] == 'Upper bound':
                                    rsr.add(forms.upper_bound(group=self.group, feature=features.solvent_access(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'Solvent access' and v[1] == 'Gaussian':
                                    rsr.add(forms.gaussian(group=self.group, feature=features.solvent_access(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                           stdev=float(v[6])))

                                if v[0] == 'Solvent access' and v[1] == 'Multiple Gaussian':
                                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.solvent_access(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                                 means=[float(v[5])],
                                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                                if v[0] == 'Solvent access' and v[1] == 'Lennard-Jones':
                                    rsr.add(forms.lennard_jones(group=self.group, feature=features.solvent_access(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]),
                                                                B=float(v[6])))

                                if v[0] == 'Solvent access' and v[1] == 'Coulomb':
                                    rsr.add(forms.coulomb(group=self.group, feature=features.solvent_access(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]),
                                                          q2=float(v[6])))

                                if v[0] == 'Solvent access' and v[1] == 'Cosine':
                                    rsr.add(forms.cosine(group=self.group, feature=features.solvent_access(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                         phase=float(v[5]),
                                                         force=float(v[6]),
                                                         period=int(v[7])))

                                # ___ density ___

                                if v[0] == 'Density' and v[1] == 'Lower bound':
                                    rsr.add(forms.lower_bound(group=self.group, feature=features.density(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'Density' and v[1] == 'Upper bound':
                                    rsr.add(forms.upper_bound(group=self.group, feature=features.density(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'Density' and v[1] == 'Gaussian':
                                    rsr.add(forms.gaussian(group=self.group, feature=features.density(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                           stdev=float(v[6])))

                                if v[0] == 'Density' and v[1] == 'Multiple Gaussian':
                                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.density(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                                 means=[float(v[5])],
                                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                                if v[0] == 'Density' and v[1] == 'Lennard-Jones':
                                    rsr.add(forms.lennard_jones(group=self.group, feature=features.density(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]),
                                                                B=float(v[6])))

                                if v[0] == 'Density' and v[1] == 'Coulomb':
                                    rsr.add(forms.coulomb(group=self.group, feature=features.density(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]),
                                                          q2=float(v[6])))

                                if v[0] == 'Density' and v[1] == 'Cosine':
                                    rsr.add(forms.cosine(group=self.group,
                                                         feature=features.density(
                                                             at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[
                                                                 8].upper()]),
                                                         phase=float(v[5]), force=float(v[6]), period=int(v[7])))

                                # ___ x coordinate ___

                                if v[0] == 'X coordinate' and v[1] == 'Lower bound':
                                    rsr.add(forms.lower_bound(group=self.group, feature=features.x_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'X coordinate' and v[1] == 'Upper bound':
                                    rsr.add(forms.upper_bound(group=self.group, feature=features.x_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'X coordinate' and v[1] == 'Gaussian':
                                    rsr.add(forms.gaussian(group=self.group, feature=features.x_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                           stdev=float(v[6])))

                                if v[0] == 'X coordinate' and v[1] == 'Multiple Gaussian':
                                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.x_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                                 means=[float(v[5])],
                                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                                if v[0] == 'X coordinate' and v[1] == 'Lennard-Jones':
                                    rsr.add(forms.lennard_jones(group=self.group, feature=features.x_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]),
                                                                B=float(v[6])))

                                if v[0] == 'X coordinate' and v[1] == 'Coulomb':
                                    rsr.add(forms.coulomb(group=self.group, feature=features.x_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]),
                                                          q2=float(v[6])))

                                if v[0] == 'X coordinate' and v[1] == 'Cosine':
                                    rsr.add(forms.cosine(group=self.group, feature=features.x_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                         phase=float(v[5]),
                                                         force=float(v[6]),
                                                         period=int(v[7])))

                                # ___ y coordinate ___

                                if v[0] == 'Y coordinate' and v[1] == 'Lower bound':
                                    rsr.add(forms.lower_bound(group=self.group, feature=features.y_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'Y coordinate' and v[1] == 'Upper bound':
                                    rsr.add(forms.upper_bound(group=self.group, feature=features.y_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'Y coordinate' and v[1] == 'Gaussian':
                                    rsr.add(forms.gaussian(group=self.group, feature=features.y_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                           stdev=float(v[6])))

                                if v[0] == 'Y coordinate' and v[1] == 'Multiple Gaussian':
                                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.y_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                                 means=[float(v[5])],
                                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                                if v[0] == 'Y coordinate' and v[1] == 'Lennard-Jones':
                                    rsr.add(forms.lennard_jones(group=self.group, feature=features.y_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]),
                                                                B=float(v[6])))

                                if v[0] == 'Y coordinate' and v[1] == 'Coulomb':
                                    rsr.add(forms.coulomb(group=self.group, feature=features.y_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]),
                                                          q2=float(v[6])))

                                if v[0] == 'Y coordinate' and v[1] == 'Cosine':
                                    rsr.add(forms.cosine(group=self.group, feature=features.y_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                         phase=float(v[5]),
                                                         force=float(v[6]),
                                                         period=int(v[7])))

                                # ___ z coordinate ___

                                if v[0] == 'Z coordinate' and v[1] == 'Lower bound':
                                    rsr.add(forms.lower_bound(group=self.group, feature=features.z_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'Z coordinate' and v[1] == 'Upper bound':
                                    rsr.add(forms.upper_bound(group=self.group, feature=features.z_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                              stdev=float(v[6])))

                                if v[0] == 'Z coordinate' and v[1] == 'Gaussian':
                                    rsr.add(forms.gaussian(group=self.group, feature=features.z_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), mean=float(v[5]),
                                                           stdev=float(v[6])))

                                if v[0] == 'Z coordinate' and v[1] == 'Multiple Gaussian':
                                    rsr.add(forms.multi_gaussian(group=self.group, feature=features.z_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                                 means=[float(v[5])],
                                                                 stdevs=[float(v[6])], weights=[float(v[7])]))

                                if v[0] == 'Z coordinate' and v[1] == 'Lennard-Jones':
                                    rsr.add(forms.lennard_jones(group=self.group, feature=features.z_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), A=float(v[5]),
                                                                B=float(v[6])))

                                if v[0] == 'Z coordinate' and v[1] == 'Coulomb':
                                    rsr.add(forms.coulomb(group=self.group, feature=features.z_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[7].upper()]), q1=float(v[5]),
                                                          q2=float(v[6])))

                                if v[0] == 'Z coordinate' and v[1] == 'Cosine':
                                    rsr.add(forms.cosine(group=self.group, feature=features.z_coordinate(
                                        at[str(v[3]).upper() + ':' + str(v[4]) + ':' + v[8].upper()]),
                                                         phase=float(v[5]),
                                                         force=float(v[6]),
                                                         period=int(v[7])))
                    # add symmetry restraints
                    elif infoFile['multi_chain'] == 'Yes' and infoFile['use_symmetry'] == 'Yes' and infoFile[
                        'symmetry'] != []:
                        # def function for residue ranges
                        for i in infoFile['symmetry']:
                            sel3 = selection(self.chains[i[1].upper()]).only_atom_types(i[0].upper())
                            sel4 = selection(self.chains[i[3].upper()]).only_atom_types(i[2].upper())
                            rsr.symmetry.append(symmetry(sel3, sel4, i[4]))
                        self.restraints.symmetry.report(3)

            a = MyModel(env, alnfile=os.path.join(self.path, self.alignFile.currentText()),
                        knowns=self.templates.currentItem().text(), sequence=self.seqId.text(), deviation=int(self.deviation.text()),
                        library_schedule=lib_schedule, assess_methods=assess_methods)
            print(self.templates.currentItem().text())
            a.starting_model = 1
            a.ending_model = int(self.numModel.text())
            a.make()
            fileL_after = [i for i in os.listdir(os.getcwd())]
            # get model files
            for i in fileL_after:
                if i not in fileL_before:
                    if os.path.exists(os.path.join(self.path, i)):
                        os.unlink(os.path.join(self.path, i))
                    shutil.move(os.path.join(os.getcwd(), i), self.path)
            self.msgLabel.setStyleSheet('color: green')
            self.msgLabel.setText('Finished')
            self.buildBut.setText('Build')
            self.buildBut.setEnabled(True)
            # display results
            self.group_model_result.show()
            # set tableView in full width
            self.table_model_result.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            #
            model_output = a.outputs
            for n in range(len(a.outputs)):
                if 'GA341 score' in a.outputs[n]:
                    model_output[n]['GA341 score'] = a.outputs[n]['GA341 score'][0]
            self.model_df = pd.DataFrame(model_output)
            self.model_df.drop(['failure', 'num', 'pdfterms'], axis=1, inplace=True)
            self.model_df = self.model_df[self.model_df.columns[::-1]]
            qt_model = pd_model.PdModel(self.model_df)
            self.table_model_result.setModel(qt_model)
        except ModuleNotFoundError:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('MODELLER not found')
            self.buildBut.setText('Build')
            self.buildBut.setEnabled(True)
        except Exception as er:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Error')
            self.buildBut.setText('Build')
            self.buildBut.setEnabled(True)
            #print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            pass

# def main():
#     app = QApplication(sys.argv)
#     form = BuildModel()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()