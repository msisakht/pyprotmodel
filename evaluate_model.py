import os
import sys
import math
import json
import collections
from Bio.SeqIO.PirIO import PirIterator
from Bio.PDB import PDBParser, PPBuilder
import matplotlib.pyplot as plt
from matplotlib import colors
import pylustrator
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QMainWindow, QApplication, QColorDialog, QHeaderView
import tpl_evaluate_model
import pd_model
import config


class EvaluateModel(QMainWindow, tpl_evaluate_model.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    pdb_fileL = []
    aln_fileL = []
    assessL = ['DOPE', 'DOPE-HR', 'Normalized DOPE', 'Normalized DOPE-HR', 'GA341']
    templateDic = {}
    addedTempColorDic = {}
    rama_types = []
    formats = plt.gcf().canvas.get_supported_filetypes_grouped()
    outlier_dic = collections.OrderedDict({'Model': [], 'Type': [], 'Chain': [], 'Residue': [], 'Residue number': []})
    countOutlier_dic = collections.OrderedDict({'Model': [], 'General': [], 'Glycine': [], 'Proline': [], 'pre-Proline': []})
    allowed_dic = collections.OrderedDict({'Model': []})

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.model.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.pdb')])
        #
        self.alignFile.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.aln')])
        self.alignFile.currentIndexChanged.connect(self.get_template)
        #
        self.assess.addItems(self.assessL)
        self.assess.currentIndexChanged.connect(self.assess_label)
        # get id
        try:
            f = open(os.path.join(self.path, self.alignFile.currentText()))
            p = list(PirIterator(f))
            ids = [i.id for i in p]
            # set id
            self.modelId.insert(0, ids[len(ids) - 1:])
        except:
            pass
        #
        self.pickColor.released.connect(self.pick_color)
        #
        self.addButton.released.connect(self.add_color)
        self.minButton.released.connect(self.min_color)
        #
        self.figSizeX.setText('12')
        self.figSizeY.setText('8')
        #
        self.lineWidth.setText('2')
        #
        self.xlabel.setText('Alignment position')
        self.ylabel.setText('DOPE per-residue score')
        #
        self.energy_save_as.addItems([k + ' (*.' + v[0] + ')' for k, v in plt.gcf().canvas.get_supported_filetypes_grouped().items()] + ['Display & edit plot'])
        self.energy_save_as.setCurrentIndex(2)
        #
        self.evalBut.released.connect(self.evaluate_energy)
        # ___ ramachandran ___
        self.rama_pdb.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        #
        self.annot.addItems(['Residue', 'Residue + Residue number', 'Residue + Residue number + Chain id', 'None'])
        self.annot.setCurrentIndex(3)
        #
        self.figSizeX_rama.setText('8')
        self.figSizeY_rama.setText('6')
        #
        self.display_outlier.addItems(['Outliers', 'Number of outliers', 'Models with outlier', 'Models without outlier'])
        self.display_outlier.currentIndexChanged.connect(self.display_outlier_format)
        #
        self.rama_save_as.addItems(['None'] + [k + ' (*.' + v[0] + ')' for k, v in plt.gcf().canvas.get_supported_filetypes_grouped().items()] + ['Display & edit plot'])
        self.rama_save_as.setCurrentIndex(0)
        #
        self.group_outlier.hide()
        #
        self.pdb_fileL = [os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')]
        self.aln_fileL = [os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.aln')]
        #
        self.evalRamaBut.released.connect(self.evaluate_rama)
        #
        self.saveCsvBut.released.connect(self.save_outlier_csv)
        self.removeBut.released.connect(self.remove_pdb_file)

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
            self.model.clear()
            self.model.addItems([os.path.basename(i) for i in os.listdir(path) if i.endswith('.pdb')])
            self.rama_pdb.clear()
            self.rama_pdb.addItems([os.path.basename(i) for i in os.listdir(path) if i.endswith('.pdb')])
            self.pdb_fileL = pdb_files
        while self.aln_fileL != aln_files:
            self.alignFile.clear()
            self.alignFile.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.aln')])
            self.aln_fileL = aln_files

    def get_template(self):
        self.templateDic = {}
        templateL = []
        self.modelId.clear()
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

                self.template_2.clear()
                self.template_2.addItems([i for i in templateL[:len(ids) - 1]])
                self.modelId.setText(ids[::-1][0])
        except:
            self.template_2.addItems([''])

    def pick_color(self):
        color = QColorDialog.getColor()
        self.tempColorEntry.setText(color.name())

    def add_color(self):
        for k, v in self.templateDic.items():
            if k == self.template_2.currentText():
                self.addedTempColorDic[k + ' (' + self.tempColorEntry.text().strip() + ')'] = v + [self.tempColorEntry.text().strip()]
        self.tempColor.clear()
        self.tempColor.addItems([k for k in self.addedTempColorDic.keys()])

    def min_color(self):
        selected = self.tempColor.selectedItems()
        for i in selected:
            del self.addedTempColorDic[i.text()]
        self.tempColor.clear()
        self.tempColor.addItems([k for k in self.addedTempColorDic.keys()])

    def assess_label(self):
        # set y label plot
        if self.assess.currentText() == 'DOPE' or self.assess.currentText() == 'DOPE-HR':
            self.ylabel.clear()
            self.ylabel.setText(self.assess.currentText() + ' per-residue score')

    @staticmethod
    def r_enumerate(seq):
        num = len(seq) - 1
        while num >= 0:
            yield num, seq[num]
            num -= 1

    @classmethod
    def get_profile(cls, profile_file, seq):
        """Read `profile_file` into a Python array, and add gaps corresponding to
           the alignment sequence `seq`."""
        # Read all non-comment and non-blank lines from the file:
        f = open(profile_file)
        vals = []
        for line in f:
            if not line.startswith('#') and len(line) > 10:
                spl = line.split()
                vals.append(float(spl[-1]))
        # Insert gaps into the profile corresponding to those in seq:
        for n, res in cls.r_enumerate(seq.residues):
            for gap in range(res.get_leading_gaps()):
                vals.insert(n, None)
        # Add a gap at position '0', so that we effectively count from 1:
        vals.insert(0, None)
        return vals

    def remove_pdb_file(self):
        try:
            for file in self.rama_pdb.selectedItems():
                os.unlink(os.path.join(self.path, file.text()))
            self.rama_pdb.clear()
            self.model.clear()
            self.rama_pdb.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
            self.model.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        except:
            pass

    def evaluate_energy(self):
        pylustrator.start()
        plt_dic = {}
        modeller_path = config.Config.get_modeller_path()
        sys.path.insert(0, modeller_path)
        try:
            from modeller import log, environ, alignment, selection
            from modeller.scripts import complete_pdb
            if self.alignFile.currentText() != 'Select':
                self.evalBut.setEnabled(False)
                self.zscore.setText('')
                self.msg.setStyleSheet('color: black')
                self.msg.setText('Processing...')
                QApplication.processEvents()
                saveName = ''
                # get models
                models = self.model.currentText()
                # get alignment file >> os.path.join(path, self.alignFile.currentText())
                # get figure size, line width & x, y labels >> from text()
                # run evaluation
                log.verbose()  # request verbose output
                env = environ()
                env.libs.topology.read(file='$(LIB)/top_heav.lib')  # read topology
                env.libs.parameters.read(file='$(LIB)/par.lib')  # read parameters

                mdl = complete_pdb(env, os.path.join(self.path, models))
                s_mdl = selection(mdl)  # all atom selection
                mdl_profile_name = os.path.join(self.path, models + '.profile')
                if self.assess.currentText() == 'DOPE' or self.assess.currentText() == 'DOPE-HR':
                    if self.assess.currentText() == 'DOPE':
                        s_mdl.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=mdl_profile_name,
                                          normalize_profile=True, smoothing_window=15)
                    elif self.assess.currentText() == 'DOPE-HR':
                        s_mdl.assess_dopehr(output='ENERGY_PROFILE NO_REPORT', file=mdl_profile_name,
                                          normalize_profile=True, smoothing_window=15)

                    # draw plot
                    if self.assess.currentText() == 'DOPE' or self.assess.currentText() == 'DOPE-HR':
                        file_format = str()
                        if self.assess.currentText() == 'DOPE':
                            saveName = models + '_dope_plot'
                        elif self.assess.currentText() == 'DOPE-HR':
                            saveName = models + '_dopehr_plot'
                        if self.energy_save_as.currentText() != 'Display & edit plot':
                            file_format = [v[0] for k, v in self.formats.items() if
                                           k + ' (*.' + v[0] + ')' == self.energy_save_as.currentText()][0]
                            self.msg.setStyleSheet('color: #000000')
                            self.msg.setText('Saving %s file...' % (file_format))
                            QApplication.processEvents()

                        a = alignment(env, file=os.path.join(self.path, self.alignFile.currentText()))
                        model = self.get_profile(mdl_profile_name, a[self.modelId.text().strip()])
                        tempProfileL = []
                        for v in self.addedTempColorDic.values():
                            # write template profile
                            tmp = complete_pdb(env, os.path.join(self.path, v[0]))
                            s_temp = selection(tmp)
                            tmp_profile_name = os.path.join(self.path, v[0] + '.profile')
                            s_temp.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=tmp_profile_name,
                                  normalize_profile=True, smoothing_window=15)
                            tempProfileL.append([tmp_profile_name, v[1], v[2]])
                        # Plot the template and model profiles in the same plot for comparison
                        plt.figure(1, figsize=(int(self.figSizeX.text().strip()), int(self.figSizeY.text().strip())))
                        plt.xlabel(self.xlabel.text().strip())
                        plt.ylabel(self.ylabel.text().strip())
                        plt.plot(model, color='red', linewidth=int(self.lineWidth.text().strip()), label='Model')
                        plt_dic['Model'] = model
                        for profile in tempProfileL:
                            template = self.get_profile(profile[0], a[profile[1]])
                            plt_dic[os.path.basename(os.path.splitext(profile[0])[0])] = template
                            plt.plot(template, color=profile[2], linewidth=int(self.lineWidth.text().strip()), label=profile[1])
                        plt.legend()
                        self.msg.setStyleSheet('color: #000000')
                        self.msg.setText('Saving plot data...')
                        QApplication.processEvents()
                        pd.DataFrame.to_csv(pd.DataFrame(dict([(k, pd.Series(v)) for k, v in plt_dic.items()])), os.path.join(self.path, saveName + '_data.csv'), index=False)
                        #
                        if self.energy_save_as.currentText() == 'Display & edit plot':
                            plt.show()
                            plt.close()
                        else:
                            plt.savefig(os.path.join(self.path, saveName + '.' + file_format))
                        print('Finished.')
                    else:
                        print('Finished.')

                elif self.assess.currentText() == 'Normalized DOPE':
                    zscore = mdl.assess_normalized_dope()
                    self.zscore.setText('Z-score: ' + str(round(zscore, 3)))
                    print('Finished.')
                elif self.assess.currentText() == 'Normalized DOPE-HR':
                    zscore = mdl.assess_normalized_dopehr()
                    self.zscore.setText('Z-score: ' + str(round(zscore, 3)))
                    print('Finished.')
                elif self.assess.currentText() == 'GA341':
                    zscore = mdl.assess_ga341()
                    self.zscore.setText('Z-score: ' + str(round(zscore, 3)))
                    print('Finished.')
                self.msg.setStyleSheet('color: green')
                self.msg.setText('Finished')
                self.evalBut.setText('Evaluate')
                self.evalBut.setEnabled(True)
            else:
                self.msg.setText('')
                QApplication.processEvents()
        except ModuleNotFoundError:
            self.msg.setStyleSheet('color: red')
            self.msg.setText('MODELLER not found')
            self.evalBut.setEnabled(True)
        except Exception as er:
            self.msg.setStyleSheet('color: red')
            self.msg.setText('Error')
            self.evalBut.setEnabled(True)
            # print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            print(er)

    # ______ ramachandran ______

    def display_outlier_format(self):
        dic = collections.OrderedDict({})
        if self.display_outlier.currentText() == 'Outliers':
            dic = self.outlier_dic
            self.outlier_count_label.setText(str(len(dic['Model'])))
        if self.display_outlier.currentText() == 'Number of outliers':
            for k, v in self.countOutlier_dic.items():
                if k == 'Model':
                    dic[k] = v
                if k in self.rama_types:
                    dic[k] = v
            self.outlier_count_label.setText('')
        if self.display_outlier.currentText() == 'Models with outlier':
            dic = self.outlier_model_dic
            self.outlier_count_label.setText(str(len(dic['Model'])))
        if self.display_outlier.currentText() == 'Models without outlier':
            dic = self.allowed_dic
            self.outlier_count_label.setText(str(len(dic['Model'])))
        self.result_df = pd.DataFrame.from_dict(dic)
        self.result_df.sort_values('Model', inplace=True)
        self.result_df.sort_index(inplace=True)
        qt_model = pd_model.PdModel(self.result_df)
        self.table_results.setModel(qt_model)
        self.table_results.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

    def save_outlier_csv(self):
        try:
            self.msglabel_csv.setStyleSheet('color: #000000')
            self.msglabel_csv.setText('Saved')
            QApplication.processEvents()
            pd.DataFrame.to_csv(self.result_df, 'Outliers.csv', index=False)
        except Exception as er:
            self.msglabel_csv.setStyleSheet('color: red')
            self.msglabel_csv.setText('Error')
            print(er)
            pass

    def evaluate_rama(self):
        pylustrator.start()
        self.rama_types = []
        rama_preferences = {}
        self.outlier_dic = collections.OrderedDict({'Model': [], 'Type': [], 'Chain': [], 'Residue': [], 'Residue number': []})
        self.countOutlier_dic = collections.OrderedDict({'Model': [], 'General': [], 'Glycine': [], 'Proline': [], 'pre-Proline': []})
        self.allowed_dic = collections.OrderedDict({'Model': []})
        self.outlier_model_dic = collections.OrderedDict({'Model': []})

        if len(self.rama_pdb.selectedItems()) > 0:
            self.msglabel_rama.setStyleSheet('color: #000000')
            self.msglabel_rama.setText('Processing...')
            QApplication.processEvents()

            try:
                if self.rama_gen.isChecked():
                    self.rama_types.append('General')
                if self.rama_gly.isChecked():
                    self.rama_types.append('Glycine')
                if self.rama_pro.isChecked():
                    self.rama_types.append('Proline')
                if self.rama_prePro.isChecked():
                    self.rama_types.append('pre-Proline')

                # General variable for the background preferences
                rama_pref_files = {
                    'General': {
                        'file': os.path.join(os.getcwd(), 'data', 'pref_general.data'),
                        'cmap': colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
                        'bounds': [0, 0.0005, 0.02, 1],
                    },
                    'Glycine': {
                        'file': os.path.join(os.getcwd(), 'data', 'pref_glycine.data'),
                        'cmap': colors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
                        'bounds': [0, 0.002, 0.02, 1],
                    },
                    'Proline': {
                        'file': os.path.join(os.getcwd(), 'data', 'pref_proline.data'),
                        'cmap': colors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
                        'bounds': [0, 0.002, 0.02, 1],
                    },
                    'pre-Proline': {
                        'file': os.path.join(os.getcwd(), 'data', 'pref_preproline.data'),
                        'cmap': colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
                        'bounds': [0, 0.002, 0.02, 1],
                    }
                }
                # get selected types
                for k, v in rama_pref_files.items():
                    if k in self.rama_types:
                        rama_preferences[k] = v
                # Read in the expected torsion angles
                __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
                rama_pref_values = {}
                for key, val in rama_preferences.items():
                    rama_pref_values[key] = np.full((360, 360), 0, dtype=np.float64)
                    with open(os.path.join(__location__, val['file'])) as fn:
                        for line in fn:
                            if not line.startswith('#'):
                                # Preference file has values for every second position only
                                rama_pref_values[key][int(float(line.split()[1])) + 180][
                                    int(float(line.split()[0])) + 180] = float(
                                    line.split()[2])
                                rama_pref_values[key][int(float(line.split()[1])) + 179][
                                    int(float(line.split()[0])) + 179] = float(
                                    line.split()[2])
                                rama_pref_values[key][int(float(line.split()[1])) + 179][
                                    int(float(line.split()[0])) + 180] = float(
                                    line.split()[2])
                                rama_pref_values[key][int(float(line.split()[1])) + 180][
                                    int(float(line.split()[0])) + 179] = float(
                                    line.split()[2])

                # Calculate the torsion angle of the inputs
                def get_plot(n, file, save=False):
                    aa_type = ''
                    plt_dic = collections.OrderedDict()
                    self.msglabel_rama.setStyleSheet('color: #000000')
                    self.msglabel_rama.setText('Processing %s/%s...' % (str(n + 1), str(len(self.rama_pdb.selectedItems()))))
                    QApplication.processEvents()
                    #
                    normals = {}
                    outliers = {}
                    for key, val in rama_preferences.items():
                        normals[key] = {'x': [], 'y': []}
                        outliers[key] = {'x': [], 'y': []}
                    #
                    structure = PDBParser().get_structure(file.text(), os.path.join(self.path, file.text()))
                    for model in structure:
                        for chain in model:
                            polypeptides = PPBuilder().build_peptides(chain)
                            for poly_index, poly in enumerate(polypeptides):
                                phi_psi = poly.get_phi_psi_list()
                                for res_index, residue in enumerate(poly):
                                    res_name = '{}'.format(residue.resname)
                                    res_num = residue.id[1]
                                    phi, psi = phi_psi[res_index]
                                    if phi and psi:
                                        if str(poly[res_index + 1].resname) == 'PRO' and 'pre-Proline' in self.rama_types:
                                            aa_type = 'pre-Proline'
                                        elif res_name == 'PRO' and 'Proline' in self.rama_types:
                                            aa_type = 'Proline'
                                        elif res_name == 'GLY' and 'Glycine' in self.rama_types:
                                            aa_type = 'Glycine'
                                        else:
                                            if 'General' in self.rama_types:
                                                aa_type = 'General'
                                        if aa_type.strip() != '' and rama_pref_values[aa_type][int(math.degrees(psi)) + 180][
                                            int(math.degrees(phi)) + 180] < rama_preferences[aa_type]['bounds'][1]:
                                            outliers[aa_type]['x'].append([math.degrees(phi), chain.id, res_name, res_num])
                                            outliers[aa_type]['y'].append([math.degrees(psi), chain.id, res_name, res_num])
                                            self.outlier_dic['Model'].append(file.text())
                                            self.outlier_dic['Type'].append(aa_type)
                                            self.outlier_dic['Chain'].append(chain.id)
                                            self.outlier_dic['Residue'].append(res_name)
                                            self.outlier_dic['Residue number'].append(res_num)
                                        else:
                                            if aa_type.strip() != '':
                                                normals[aa_type]['x'].append([math.degrees(phi), chain.id, res_name, res_num])
                                                normals[aa_type]['y'].append([math.degrees(psi), chain.id, res_name, res_num])

                    # get model with & without outlier
                    if file.text() not in self.outlier_dic['Model'] and file.text() not in self.allowed_dic['Model']:
                        self.allowed_dic['Model'].append(file.text())
                    if file.text() in self.outlier_dic['Model'] and file.text() not in self.outlier_model_dic['Model']:
                        self.outlier_model_dic['Model'].append(file.text())
                    # count outliers
                    self.countOutlier_dic['Model'].append(file.text())
                    if 'General' in self.rama_types:
                        self.countOutlier_dic['General'].append(len(outliers['General']['x']))
                    if 'Glycine' in self.rama_types:
                        self.countOutlier_dic['Glycine'].append(len(outliers['Glycine']['x']))
                    if 'Proline' in self.rama_types:
                        self.countOutlier_dic['Proline'].append(len(outliers['Proline']['x']))
                    if 'pre-Proline' in self.rama_types:
                        self.countOutlier_dic['pre-Proline'].append(len(outliers['pre-Proline']['x']))
                    # Generate the plots
                    if self.rama_save_as.currentText() != 'None':
                        plt.figure(figsize=(int(self.figSizeX_rama.text()), int(self.figSizeY_rama.text())))
                        for idx, (key, val) in enumerate(sorted(rama_preferences.items(), key=lambda x: x[0].lower())):
                            annotatL = []
                            if len(rama_preferences) == 4:
                                plt.subplot(2, 2, idx + 1)
                            if len(rama_preferences) == 3:
                                plt.subplot(1, 3, idx + 1)
                            if len(rama_preferences) == 2:
                                plt.subplot(1, 2, idx + 1)
                            if len(rama_preferences) == 1:
                                plt.subplot(1, 1, 1)
                            plt.title(key)
                            plt.imshow(rama_pref_values[key], cmap=rama_preferences[key]['cmap'],
                                       norm=colors.BoundaryNorm(rama_preferences[key]['bounds'], rama_preferences[key]['cmap'].N),
                                       extent=(-180, 180, 180, -180))
                            plt.scatter([i[0] for i in normals[key]['x']], [i[0] for i in normals[key]['y']])
                            plt.scatter([i[0] for i in outliers[key]['x']], [i[0] for i in outliers[key]['y']], color='red')
                            plt_dic[key + '_normal_x'] = [i[0] for i in normals[key]['x']]
                            plt_dic[key + '_normal_y'] = [i[0] for i in normals[key]['y']]
                            plt_dic[key + '_outlier_x'] = [i[0] for i in outliers[key]['x']]
                            plt_dic[key + '_outlier_y'] = [i[0] for i in outliers[key]['y']]
                            # plt_dic['preference_values'] = rama_pref_values[key]
                            # add annotation
                            for m, txt in enumerate(outliers[key]['x']):
                                if self.annot.currentText() == 'Residue':
                                    annotat = [i[2] for i in outliers[key]['x']][m]
                                elif self.annot.currentText() == 'Residue + Residue number':
                                    annotat = [i[2] + str(i[3]) for i in outliers[key]['x']][m]
                                elif self.annot.currentText() == 'Residue + Residue number + Chain id':
                                    annotat = [i[2] + str(i[3]) + '(' + i[1] + ')' for i in outliers[key]['x']][m]
                                else:
                                    annotat = ''
                                plt.annotate(annotat, ([i[0] for i in outliers[key]['x']][m] + 2., [i[0] for i in outliers[key]['y']][m] + 2.))
                                annotatL.append(annotat)
                            plt_dic[key + '_outlier_annotation'] = annotatL
                            plt.xlim([-180, 180])
                            plt.ylim([-180, 180])
                            plt.xticks([-180, 180])
                            plt.yticks([-180, 180])
                            plt.plot([-180, 180], [0, 0], color='black')
                            plt.plot([0, 0], [-180, 180], color='black')
                            plt.locator_params(axis='x', nbins=7)
                            plt.xlabel(r'$\phi$')
                            plt.ylabel(r'$\psi$')
                            plt.grid()
                        plt.tight_layout()
                        # save plot to file
                        if save:
                            file_format = [v[0] for k, v in self.formats.items() if k + ' (*.' + v[0] + ')' == self.rama_save_as.currentText()][0]
                            if self.rama_gen.isChecked():
                                rama_gen = 'general_'
                            else:
                                rama_gen = ''
                            if self.rama_gly.isChecked():
                                rama_gly = 'gly_'
                            else:
                                rama_gly = ''
                            if self.rama_pro.isChecked():
                                rama_pro = 'pro_'
                            else:
                                rama_pro = ''
                            if self.rama_prePro.isChecked():
                                rama_prePro = 'prePro_'
                            else:
                                rama_prePro = ''
                            file_name = rama_gen + rama_gly + rama_pro + rama_prePro
                            plt.savefig(os.path.join(self.path, os.path.splitext(file.text())[0] + '_ramachandran_' + file_name + '.' + file_format))
                        else:
                            plt.show()
                            plt.close()
                        # save plot data
                        pd.DataFrame.to_csv(pd.DataFrame(dict([(k, pd.Series(v)) for k, v in plt_dic.items()])),
                                            os.path.join(self.path, os.path.splitext(file.text())[0] +
                                                         '_PhiPsi_plot_data.csv'), index=False)
                    # show outliers
                    if len(self.outlier_dic) > 0 or len(self.allowed_dic) > 0:
                        self.group_outlier.show()
                        dic = collections.OrderedDict()
                        if self.display_outlier.currentText() == 'Outliers':
                            dic = self.outlier_dic
                            self.outlier_count_label.setText(str(len(dic['Model'])))
                        if self.display_outlier.currentText() == 'Number of outliers':
                            for k, v in self.countOutlier_dic.items():
                                if k == 'Model':
                                    dic[k] = v
                                if k in self.rama_types:
                                    dic[k] = v
                            self.outlier_count_label.setText('')
                        if self.display_outlier.currentText() == 'Models with outlier':
                            dic = self.outlier_model_dic
                            self.outlier_count_label.setText(str(len(dic['Model'])))
                        if self.display_outlier.currentText() == 'Models without outlier':
                            dic = self.allowed_dic
                            self.outlier_count_label.setText(str(len(dic['Model'])))
                        self.result_df = pd.DataFrame.from_dict(dic)
                        qt_model = pd_model.PdModel(self.result_df)
                        self.table_results.setModel(qt_model)
                        self.table_results.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
                    self.msglabel_rama.setStyleSheet('color: green')
                    self.msglabel_rama.setText('Finished')

                if self.rama_save_as.currentText() == 'Display & edit plot':
                    get_plot(0, self.rama_pdb.selectedItems()[0])
                else:
                    for n, file in enumerate(self.rama_pdb.selectedItems()):
                        get_plot(n, file, True)
            except Exception as er:
                print(er)
                self.msglabel_rama.setStyleSheet('color: red')
                self.msglabel_rama.setText('Error')
                # print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))


# def main():
#     app = QApplication(sys.argv)
#     form = EvaluateModel()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()