import os
import sys
import math
import json
import shutil
import threading
import collections
from Bio.SeqIO.PirIO import PirIterator
from Bio.PDB import PDBParser, PPBuilder
import matplotlib.pyplot as plt
from matplotlib import colors
import pylustrator
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QMainWindow, QApplication, QColorDialog, QHeaderView, QFileDialog
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
        self.browseModel1.released.connect(self.browse_model_th)
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
        self.browsModel2.released.connect(self.browse_model_th)
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

    def browse_model_th(self):
        threading.Thread(target=self.browse_model).start()

    def browse_model(self):
        try:
            user_pdb = QFileDialog.getOpenFileNames()[0]
            if user_pdb:
                for file in user_pdb:
                    shutil.copy(file, self.path)
                self.model.clear()
                self.rama_pdb.clear()
                self.rama_pdb.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
                self.model.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
                self.model.setCurrentText(os.path.basename(user_pdb[0]))
        except:
            pass

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
        dic = collections.OrderedDict()

        if self.display_outlier.currentText() == 'Outliers':
            dic = self.outlier_dic
            self.outlier_count_label.setText(str(len(dic['Model'])))

        elif self.display_outlier.currentText() == 'Number of outliers':
            for k, v in self.countOutlier_dic.items():
                if k == 'Model' or k in self.rama_types:
                    dic[k] = v
            self.outlier_count_label.setText('')

        elif self.display_outlier.currentText() == 'Models with outlier':
            dic = self.outlier_model_dic
            self.outlier_count_label.setText(str(len(dic['Model'])))

        elif self.display_outlier.currentText() == 'Models without outlier':
            dic = self.allowed_dic
            self.outlier_count_label.setText(str(len(dic['Model'])))

        self.result_df = pd.DataFrame.from_dict(dic)

        if not self.result_df.empty and 'Model' in self.result_df.columns:
            self.result_df.sort_values('Model', inplace=True)
            self.result_df.reset_index(drop=True, inplace=True)

        qt_model = pd_model.PdModel(self.result_df)
        self.table_results.setModel(qt_model)
        self.table_results.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)


    def save_outlier_csv(self):
        try:
            if not hasattr(self, 'result_df') or self.result_df.empty:
                self.msglabel_csv.setStyleSheet('color: red')
                self.msglabel_csv.setText('No data to save')
                return

            self.msglabel_csv.setStyleSheet('color: #000000')
            self.msglabel_csv.setText('Saving...')
            QApplication.processEvents()

            self.result_df.to_csv(
                os.path.join(self.path, 'Outliers.csv'),
                index=False
            )

            self.msglabel_csv.setStyleSheet('color: green')
            self.msglabel_csv.setText('Saved')

        except Exception as er:
            self.msglabel_csv.setStyleSheet('color: red')
            self.msglabel_csv.setText('Error')
            print(er)


    def evaluate_rama(self):
        from matplotlib.lines import Line2D

        pylustrator.start()

        self.rama_types = []
        self.outlier_dic = collections.OrderedDict({
            'Model': [],
            'Type': [],
            'Chain': [],
            'Residue': [],
            'Residue number': []
        })
        self.countOutlier_dic = collections.OrderedDict({
            'Model': [],
            'General': [],
            'Glycine': [],
            'Proline': [],
            'pre-Proline': []
        })
        self.allowed_dic = collections.OrderedDict({'Model': []})
        self.outlier_model_dic = collections.OrderedDict({'Model': []})

        selected_pdb_files = self.rama_pdb.selectedItems()

        if len(selected_pdb_files) == 0:
            self.msglabel_rama.setStyleSheet('color: red')
            self.msglabel_rama.setText('Select at least one PDB model')
            return

        if self.rama_gen.isChecked():
            self.rama_types.append('General')

        if self.rama_gly.isChecked():
            self.rama_types.append('Glycine')

        if self.rama_pro.isChecked():
            self.rama_types.append('Proline')

        if self.rama_prePro.isChecked():
            self.rama_types.append('pre-Proline')

        if len(self.rama_types) == 0:
            self.msglabel_rama.setStyleSheet('color: red')
            self.msglabel_rama.setText('Select at least one residue category')
            return

        self.msglabel_rama.setStyleSheet('color: #000000')
        self.msglabel_rama.setText('Processing...')
        QApplication.processEvents()

        try:
            module_dir = os.path.dirname(os.path.abspath(__file__))

            # Background map settings:
            # White = disallowed region
            # Light color = allowed region
            # Dark color = favoured region
            rama_pref_files = {
                'General': {
                    'file': os.path.join(module_dir, 'data', 'pref_general.data'),
                    'cmap': colors.ListedColormap([
                        '#FFFFFF',
                        '#CBE8F6',
                        '#67B7DD'
                    ]),
                    'bounds': [0, 0.0005, 0.02, 1]
                },
                'Glycine': {
                    'file': os.path.join(module_dir, 'data', 'pref_glycine.data'),
                    'cmap': colors.ListedColormap([
                        '#FFFFFF',
                        '#FCE1B0',
                        '#F3AB4F'
                    ]),
                    'bounds': [0, 0.002, 0.02, 1]
                },
                'Proline': {
                    'file': os.path.join(module_dir, 'data', 'pref_proline.data'),
                    'cmap': colors.ListedColormap([
                        '#FFFFFF',
                        '#D8EDC9',
                        '#72B85B'
                    ]),
                    'bounds': [0, 0.002, 0.02, 1]
                },
                'pre-Proline': {
                    'file': os.path.join(module_dir, 'data', 'pref_preproline.data'),
                    'cmap': colors.ListedColormap([
                        '#FFFFFF',
                        '#D6EAF8',
                        '#6EAED4'
                    ]),
                    'bounds': [0, 0.002, 0.02, 1]
                }
            }

            rama_preferences = {
                key: value
                for key, value in rama_pref_files.items()
                if key in self.rama_types
            }

            # Read Ramachandran preference maps.
            rama_pref_values = {}

            for key, val in rama_preferences.items():
                rama_pref_values[key] = np.zeros((360, 360), dtype=np.float64)

                with open(val['file'], 'r') as fn:
                    for line in fn:
                        if line.startswith('#'):
                            continue

                        parts = line.split()

                        if len(parts) < 3:
                            continue

                        phi_value = int(float(parts[0]))
                        psi_value = int(float(parts[1]))
                        preference_value = float(parts[2])

                        phi_index = phi_value + 180
                        psi_index = psi_value + 180

                        valid_indices = [
                            (psi_index, phi_index),
                            (psi_index - 1, phi_index - 1),
                            (psi_index - 1, phi_index),
                            (psi_index, phi_index - 1)
                        ]

                        for row, col in valid_indices:
                            if 0 <= row < 360 and 0 <= col < 360:
                                rama_pref_values[key][row, col] = preference_value

            def angle_to_preference_index(angle_in_radians):
                """
                Convert a torsion angle in radians to a valid preference-map index.
                The map stores values over -180 to +179 degrees.
                """
                angle_in_degrees = int(math.degrees(angle_in_radians))
                angle_in_degrees = int(np.clip(angle_in_degrees, -180, 179))
                return angle_in_degrees + 180

            def get_plot(model_number, pdb_item, save=False):
                plt_dic = collections.OrderedDict()

                self.msglabel_rama.setStyleSheet('color: #000000')
                self.msglabel_rama.setText(
                    'Processing %s/%s...' % (
                        str(model_number + 1),
                        str(len(selected_pdb_files))
                    )
                )
                QApplication.processEvents()

                normals = {}
                outliers = {}

                for key in rama_preferences:
                    normals[key] = {'x': [], 'y': []}
                    outliers[key] = {'x': [], 'y': []}

                pdb_name = pdb_item.text()
                pdb_path = os.path.join(self.path, pdb_name)

                structure = PDBParser(QUIET=True).get_structure(
                    pdb_name,
                    pdb_path
                )

                for model in structure:
                    for chain in model:
                        polypeptides = PPBuilder().build_peptides(chain)

                        for poly in polypeptides:
                            phi_psi_list = poly.get_phi_psi_list()

                            for residue_index, residue in enumerate(poly):
                                phi, psi = phi_psi_list[residue_index]

                                # Do not ignore valid 0° torsion angles.
                                if phi is None or psi is None:
                                    continue

                                residue_name = residue.resname
                                residue_number = residue.id[1]

                                # Reset residue type at every residue.
                                aa_type = None

                                next_residue_is_proline = (
                                    residue_index + 1 < len(poly)
                                    and poly[residue_index + 1].resname == 'PRO'
                                )

                                if next_residue_is_proline and 'pre-Proline' in self.rama_types:
                                    aa_type = 'pre-Proline'

                                elif residue_name == 'PRO' and 'Proline' in self.rama_types:
                                    aa_type = 'Proline'

                                elif residue_name == 'GLY' and 'Glycine' in self.rama_types:
                                    aa_type = 'Glycine'

                                elif 'General' in self.rama_types:
                                    aa_type = 'General'

                                if aa_type is None:
                                    continue

                                phi_degree = math.degrees(phi)
                                psi_degree = math.degrees(psi)

                                phi_index = angle_to_preference_index(phi)
                                psi_index = angle_to_preference_index(psi)

                                preference_value = rama_pref_values[aa_type][
                                    psi_index,
                                    phi_index
                                ]

                                is_outlier = (
                                    preference_value
                                    < rama_preferences[aa_type]['bounds'][1]
                                )

                                residue_record_x = [
                                    phi_degree,
                                    chain.id,
                                    residue_name,
                                    residue_number
                                ]
                                residue_record_y = [
                                    psi_degree,
                                    chain.id,
                                    residue_name,
                                    residue_number
                                ]

                                if is_outlier:
                                    outliers[aa_type]['x'].append(residue_record_x)
                                    outliers[aa_type]['y'].append(residue_record_y)

                                    self.outlier_dic['Model'].append(pdb_name)
                                    self.outlier_dic['Type'].append(aa_type)
                                    self.outlier_dic['Chain'].append(chain.id)
                                    self.outlier_dic['Residue'].append(residue_name)
                                    self.outlier_dic['Residue number'].append(
                                        residue_number
                                    )

                                else:
                                    normals[aa_type]['x'].append(residue_record_x)
                                    normals[aa_type]['y'].append(residue_record_y)

                # Classify models with versus without outliers.
                if (
                    pdb_name not in self.outlier_dic['Model']
                    and pdb_name not in self.allowed_dic['Model']
                ):
                    self.allowed_dic['Model'].append(pdb_name)

                if (
                    pdb_name in self.outlier_dic['Model']
                    and pdb_name not in self.outlier_model_dic['Model']
                ):
                    self.outlier_model_dic['Model'].append(pdb_name)

                # Count outliers for each category.
                self.countOutlier_dic['Model'].append(pdb_name)

                for rama_type in self.rama_types:
                    self.countOutlier_dic[rama_type].append(
                        len(outliers[rama_type]['x'])
                    )

                # Generate and save the Ramachandran figure.
                if self.rama_save_as.currentText() != 'None':
                    ordered_preferences = sorted(
                        rama_preferences.items(),
                        key=lambda item: item[0].casefold()
                    )

                    number_of_panels = len(ordered_preferences)

                    if number_of_panels == 1:
                        number_of_rows = 1
                        number_of_columns = 1

                    elif number_of_panels == 2:
                        number_of_rows = 1
                        number_of_columns = 2

                    elif number_of_panels == 3:
                        number_of_rows = 1
                        number_of_columns = 3

                    else:
                        number_of_rows = 2
                        number_of_columns = 2

                    try:
                        requested_width = float(
                            self.figSizeX_rama.text().strip()
                        )
                    except ValueError:
                        requested_width = 8.0

                    try:
                        requested_height = float(
                            self.figSizeY_rama.text().strip()
                        )
                    except ValueError:
                        requested_height = 6.0

                    # Ensure sufficient space for larger, publication-readable text.
                    figure_width = max(
                        requested_width,
                        number_of_columns * 5.8
                    )
                    figure_height = max(
                        requested_height,
                        number_of_rows * 5.8 + 0.6
                    )

                    with plt.rc_context({
                        'font.family': 'sans-serif',
                        'font.sans-serif': ['Arial', 'DejaVu Sans'],
                        'font.size': 14,
                        'axes.labelcolor': '#111111',
                        'axes.titlecolor': '#111111',
                        'xtick.color': '#111111',
                        'ytick.color': '#111111'
                    }):
                        fig, axes = plt.subplots(
                            number_of_rows,
                            number_of_columns,
                            figsize=(figure_width, figure_height),
                            squeeze=False
                        )

                        axes_flat = axes.flatten()

                        # Major tick positions requested for both axes.
                        torsion_ticks = np.arange(-180, 181, 60)

                        for panel_index, (key, val) in enumerate(
                            ordered_preferences
                        ):
                            ax = axes_flat[panel_index]

                            normal_x = [
                                entry[0]
                                for entry in normals[key]['x']
                            ]
                            normal_y = [
                                entry[0]
                                for entry in normals[key]['y']
                            ]
                            outlier_x = [
                                entry[0]
                                for entry in outliers[key]['x']
                            ]
                            outlier_y = [
                                entry[0]
                                for entry in outliers[key]['y']
                            ]

                            # Plot preference map.
                            ax.imshow(
                                rama_pref_values[key],
                                cmap=val['cmap'],
                                norm=colors.BoundaryNorm(
                                    val['bounds'],
                                    val['cmap'].N
                                ),
                                extent=(-180, 180, -180, 180),
                                origin='lower',
                                interpolation='nearest',
                                aspect='equal',
                                zorder=0
                            )

                            # Observed residues within allowed/favoured regions.
                            if len(normal_x) > 0:
                                ax.scatter(
                                    normal_x,
                                    normal_y,
                                    s=54,
                                    c='#1F4E79',
                                    edgecolors='#FFFFFF',
                                    linewidths=0.65,
                                    alpha=0.95,
                                    zorder=3
                                )

                            # Outlier residues.
                            if len(outlier_x) > 0:
                                ax.scatter(
                                    outlier_x,
                                    outlier_y,
                                    s=60,
                                    c='#C62828',
                                    edgecolors='#FFFFFF',
                                    linewidths=0.65,
                                    alpha=1.0,
                                    zorder=4
                                )

                            # Save plot data for this residue category.
                            plt_dic[key + '_normal_x'] = normal_x
                            plt_dic[key + '_normal_y'] = normal_y
                            plt_dic[key + '_outlier_x'] = outlier_x
                            plt_dic[key + '_outlier_y'] = outlier_y

                            annotation_list = []

                            for point_index, residue_data in enumerate(
                                outliers[key]['x']
                            ):
                                if self.annot.currentText() == 'Residue':
                                    annotation = residue_data[2]

                                elif (
                                    self.annot.currentText()
                                    == 'Residue + Residue number'
                                ):
                                    annotation = (
                                        residue_data[2]
                                        + str(residue_data[3])
                                    )

                                elif (
                                    self.annot.currentText()
                                    == 'Residue + Residue number + Chain id'
                                ):
                                    annotation = (
                                        residue_data[2]
                                        + str(residue_data[3])
                                        + ' ('
                                        + residue_data[1]
                                        + ')'
                                    )

                                else:
                                    annotation = ''

                                if annotation:
                                    ax.annotate(
                                        annotation,
                                        xy=(
                                            outlier_x[point_index],
                                            outlier_y[point_index]
                                        ),
                                        xytext=(5, 5),
                                        textcoords='offset points',
                                        fontsize=11,
                                        fontweight='semibold',
                                        color='#7A1F1F',
                                        zorder=5
                                    )

                                annotation_list.append(annotation)

                            plt_dic[
                                key + '_outlier_annotation'
                            ] = annotation_list

                            # Titles, axis labels, units, ticks, and visual style.
                            ax.set_title(
                                key,
                                fontsize=18,
                                fontweight='semibold',
                                pad=10
                            )

                            ax.set_xlim(-180, 180)
                            ax.set_ylim(-180, 180)

                            ax.set_xticks(torsion_ticks)
                            ax.set_yticks(torsion_ticks)

                            ax.set_xlabel(
                                r'$\phi$ torsion angle (°)',
                                fontsize=15,
                                fontweight='semibold',
                                labelpad=8
                            )

                            ax.set_ylabel(
                                r'$\psi$ torsion angle (°)',
                                fontsize=15,
                                fontweight='semibold',
                                labelpad=8
                            )

                            ax.tick_params(
                                axis='both',
                                which='major',
                                labelsize=13,
                                width=1.2,
                                length=5
                            )

                            # Zero-angle reference axes; these are not grid lines.
                            ax.axhline(
                                0,
                                color='#222222',
                                linewidth=1.15,
                                zorder=2
                            )
                            ax.axvline(
                                0,
                                color='#222222',
                                linewidth=1.15,
                                zorder=2
                            )

                            # Explicitly remove all grid lines.
                            ax.grid(False)

                            for spine in ax.spines.values():
                                spine.set_linewidth(1.25)
                                spine.set_color('#222222')

                        # Hide unused subplot axes, if any.
                        for unused_axis in axes_flat[number_of_panels:]:
                            unused_axis.set_visible(False)

                        # Shared legend for observed residues and outliers.
                        legend_handles = [
                            Line2D(
                                [0],
                                [0],
                                marker='o',
                                color='w',
                                label='Residues within allowed regions',
                                markerfacecolor='#1F4E79',
                                markeredgecolor='#FFFFFF',
                                markeredgewidth=0.65,
                                markersize=9
                            ),
                            Line2D(
                                [0],
                                [0],
                                marker='o',
                                color='w',
                                label='Outlier residues',
                                markerfacecolor='#C62828',
                                markeredgecolor='#FFFFFF',
                                markeredgewidth=0.65,
                                markersize=9
                            )
                        ]

                        fig.legend(
                            handles=legend_handles,
                            loc='lower center',
                            ncol=2,
                            frameon=False,
                            fontsize=12,
                            bbox_to_anchor=(0.5, 0.01)
                        )

                        fig.tight_layout(rect=(0, 0.075, 1, 1))

                        if save:
                            file_format = [
                                value[0]
                                for name, value in self.formats.items()
                                if (
                                    name + ' (*.' + value[0] + ')'
                                    == self.rama_save_as.currentText()
                                )
                            ][0]

                            rama_gen = (
                                'general_'
                                if self.rama_gen.isChecked()
                                else ''
                            )
                            rama_gly = (
                                'gly_'
                                if self.rama_gly.isChecked()
                                else ''
                            )
                            rama_pro = (
                                'pro_'
                                if self.rama_pro.isChecked()
                                else ''
                            )
                            rama_prepro = (
                                'prePro_'
                                if self.rama_prePro.isChecked()
                                else ''
                            )

                            plot_suffix = (
                                rama_gen
                                + rama_gly
                                + rama_pro
                                + rama_prepro
                            )

                            figure_path = os.path.join(
                                self.path,
                                os.path.splitext(pdb_name)[0]
                                + '_ramachandran_'
                                + plot_suffix
                                + '.'
                                + file_format
                            )

                            fig.savefig(
                                figure_path,
                                dpi=300,
                                bbox_inches='tight',
                                facecolor='white'
                            )

                            plt.close(fig)

                        else:
                            plt.show()
                            plt.close(fig)

                    # Save numerical Ramachandran data.
                    pd.DataFrame(
                        dict(
                            [
                                (key, pd.Series(value))
                                for key, value in plt_dic.items()
                            ]
                        )
                    ).to_csv(
                        os.path.join(
                            self.path,
                            os.path.splitext(pdb_name)[0]
                            + '_PhiPsi_plot_data.csv'
                        ),
                        index=False
                    )

                # Display outlier information in the table.
                if (
                    len(self.outlier_dic['Model']) > 0
                    or len(self.allowed_dic['Model']) > 0
                ):
                    self.group_outlier.show()
                    self.display_outlier_format()

                self.msglabel_rama.setStyleSheet('color: green')
                self.msglabel_rama.setText('Finished')

            if self.rama_save_as.currentText() == 'Display & edit plot':
                get_plot(0, selected_pdb_files[0], save=False)

            else:
                for model_number, pdb_item in enumerate(selected_pdb_files):
                    get_plot(model_number, pdb_item, save=True)

        except Exception as er:
            print(er)
            self.msglabel_rama.setStyleSheet('color: red')
            self.msglabel_rama.setText('Error')


# def main():
#     app = QApplication(sys.argv)
#     form = EvaluateModel()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()