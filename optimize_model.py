import os
import sys
import json
import threading
from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog
import tpl_optimize_model


class OptimizeModel(QMainWindow, tpl_optimize_model.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    pdb_fileL = []
    rsr_fileL = []
    optimize_methodL = ['Conjugate gradients', 'Quasi newton', 'Molecular dynamics', '1 & 3']
    YesNoL = ['Yes', 'No']
    physicalFeatL = ['Bond', 'Angle', 'Dihedral', 'Improper', 'Soft_Sphere', 'Lennard_Jones', 'Coulomb', 'H_Bond',
                     'CA_Distance', 'N_O_Distance', 'Phi_Dihedral', 'Psi_Dihedral', 'Omega_Dihedral', 'Chi1_Dihedral',
                     'Chi2_Dihedral', 'Chi3_Dihedral', 'Chi4_Dihedral', 'Disulfide_Distance', 'Disulfide_Angle',
                     'Disulfide_Dihedral', 'Lower_Distance', 'Upper_Distance', 'SD_MN_Distance', 'Chi5_Dihedral',
                     'Phi_Psi_Dihedral', 'SD_SD_Distance', 'XY_Distance', 'NMR_Distance', 'NMR_Distance2',
                     'Min_Distance', 'Nonbond_Spline', 'Accessibility', 'Density', 'Absposition', 'Dihedral_Diff',
                     'SAXS', 'Symmetry']
    get_physicalL = []
    #
    optimize_SC_methodL = ['SCWRL4']
    frame_file = ''
    seq_file = ''
    param_file = ''

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.optimizeMethod.addItems(self.optimize_methodL)
        self.optimizeMethod.currentTextChanged.connect(self.check_method)
        #
        self.model.addItems([i for i in os.listdir(self.path) if i.endswith('.pdb')])
        self.model.setCurrentRow(0)
        #
        self.default_2.setText('1.0')
        self.bond.setText('0')
        self.angle.setText('0')
        self.dihedral.setText('0')
        self.improper.setText('0')
        self.soft_sphere.setText('0')
        self.lennard_jones.setText('0')
        self.coulomb.setText('0')
        self.h_bond.setText('0')
        self.ca_distance.setText('0')
        self.n_o_distance.setText('0')
        self.phi_dihedral.setText('0')
        self.psi_dihedral.setText('0')
        self.omega_dihedral.setText('0')
        self.chi1_dihedral.setText('0')
        self.chi2_dihedral.setText('0')
        self.chi3_dihedral.setText('0')
        self.chi4_dihedral.setText('0')
        self.chi5_dihedral.setText('0')
        self.disulfide_distance.setText('0')
        self.disulfide_angle.setText('0')
        self.disulfide_dihedral.setText('0')
        self.lower_distance.setText('0')
        self.upper_distance.setText('0')
        self.sd_mn_distance.setText('0')
        self.phi_psi_dihedral.setText('0')
        self.sd_sd_distance.setText('0')
        self.xy_distance.setText('0')
        self.nmr_distance.setText('0')
        self.nmr_distance2.setText('0')
        self.min_distance.setText('0')
        self.nonbond_spline.setText('0')
        self.accessibility.setText('0')
        self.density.setText('0')
        self.absposition.setText('0')
        self.dihedral_diff.setText('0')
        self.SAXS.setText('0')
        self.symmetry.setText('0')
        #
        self.restraint.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.rsr')])
        #
        self.minAtomShift.setText('0.01')
        #
        self.maxAtomShift.setText('100.0')
        #
        self.capAtomShift.setText('0.2')
        #
        self.mdTimeStep.setText('4.0')
        #
        self.temperature.setText('293.0')
        #
        self.maxIter.setText('20')
        #
        self.From.setText('0')
        self.To.setText('99999')
        #
        self.optBackBoneBut.released.connect(self.optimize_backBone_th)
        self.backBone_remove.released.connect(self.remove_backBone_pdb)

    def enterEvent(self, QEvent):
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']
        self.path = path

        pdb_files = []
        rsr_files = []
        for i in os.listdir(self.path):
            if i.endswith('.pdb') and i not in pdb_files:
                pdb_files.append(i)
            if i.endswith('.rsr') and i not in rsr_files:
                rsr_files.append(i)
        while self.pdb_fileL != pdb_files:
            self.model.clear()
            self.model.addItems([i for i in os.listdir(path) if i.endswith('.pdb')])
            self.pdb_fileL = pdb_files
        while self.rsr_fileL != rsr_files:
            self.restraint.clear()
            self.restraint.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.rsr')])
            self.rsr_fileL = rsr_files

    def check_method(self):
        if self.optimizeMethod.currentText() == 'Conjugate gradients':
            self.maxAtomShift.setEnabled(False)
            self.maxAtomShiftLabel.setEnabled(False)
            self.capAtomShift.setEnabled(False)
            self.capAtomShiftLabel.setEnabled(False)
            self.mdTimeStep.setEnabled(False)
            self.mdTimeStepLabel.setEnabled(False)
            self.temperature.setEnabled(False)
            self.temperatureLabel.setEnabled(False)
        if self.optimizeMethod.currentText() == 'Quasi newton':
            self.maxAtomShift.setEnabled(True)
            self.maxAtomShiftLabel.setEnabled(True)
            self.capAtomShift.setEnabled(False)
            self.capAtomShiftLabel.setEnabled(False)
            self.mdTimeStep.setEnabled(False)
            self.mdTimeStepLabel.setEnabled(False)
            self.temperature.setEnabled(False)
            self.temperatureLabel.setEnabled(False)
        if self.optimizeMethod.currentText() == 'Molecular dynamics':
            self.maxAtomShift.setEnabled(False)
            self.maxAtomShiftLabel.setEnabled(False)
            self.capAtomShift.setEnabled(True)
            self.capAtomShiftLabel.setEnabled(True)
            self.mdTimeStep.setEnabled(True)
            self.mdTimeStepLabel.setEnabled(True)
            self.temperature.setEnabled(True)
            self.temperatureLabel.setEnabled(True)
        if self.optimizeMethod.currentText() == '1 & 3':
            self.maxAtomShift.setEnabled(False)
            self.maxAtomShiftLabel.setEnabled(False)
            self.capAtomShift.setEnabled(True)
            self.capAtomShiftLabel.setEnabled(True)
            self.mdTimeStep.setEnabled(True)
            self.mdTimeStepLabel.setEnabled(True)
            self.temperature.setEnabled(True)
            self.temperatureLabel.setEnabled(True)

    def remove_backBone_pdb(self):
        try:
            for file in self.model.selectedItems():
                os.unlink(os.path.join(self.path, file.text()))
            self.model.clear()
            self.model.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
        except:
            pass

    def optimize_backBone_th(self):
        threading.Thread(target=self.optimize_backBone).start()
        self.optBackBoneBut.setEnabled(False)

    def optimize_backBone(self):
        self.msg_1.setText('Processing...')
        modeller_path = self.config.Config.get_modeller_path()
        sys.path.insert(0, modeller_path)
        try:
            from modeller import environ, selection, physical
            from modeller.scripts import complete_pdb
            from modeller.optimizers import conjugate_gradients, quasi_newton, molecular_dynamics, actions
            if len(self.model.currentItem().text()) > 0:
                # get physical features
                bond = float(self.default_2.text().strip())
                angle = float(self.default_2.text().strip())
                dihedral = float(self.default_2.text().strip())
                improper = float(self.default_2.text().strip())
                soft_sphere = float(self.default_2.text().strip())
                lennard_jones = float(self.default_2.text().strip())
                coulomb = float(self.default_2.text().strip())
                h_bond = float(self.default_2.text().strip())
                ca_distance = float(self.default_2.text().strip())
                n_o_distance = float(self.default_2.text().strip())
                phi_dihedral = float(self.default_2.text().strip())
                psi_dihedral = float(self.default_2.text().strip())
                omega_dihedral = float(self.default_2.text().strip())
                chi1_dihedral = float(self.default_2.text().strip())
                chi2_dihedral = float(self.default_2.text().strip())
                chi3_dihedral = float(self.default_2.text().strip())
                chi4_dihedral = float(self.default_2.text().strip())
                disulfide_distance = float(self.default_2.text().strip())
                disulfide_angle = float(self.default_2.text().strip())
                disulfide_dihedral = float(self.default_2.text().strip())
                lower_distance = float(self.default_2.text().strip())
                upper_distance = float(self.default_2.text().strip())
                sd_mn_distance = float(self.default_2.text().strip())
                chi5_dihedral = float(self.default_2.text().strip())
                phi_psi_dihedral = float(self.default_2.text().strip())
                sd_sd_distance = float(self.default_2.text().strip())
                xy_distance = float(self.default_2.text().strip())
                nmr_distance = float(self.default_2.text().strip())
                nmr_distance2 = float(self.default_2.text().strip())
                min_distance = float(self.default_2.text().strip())
                nonbond_spline = float(self.default_2.text().strip())
                accessibility = float(self.default_2.text().strip())
                density = float(self.default_2.text().strip())
                absposition = float(self.default_2.text().strip())
                dihedral_diff = float(self.default_2.text().strip())
                saxs = float(self.default_2.text().strip())
                symmetry = float(self.default_2.text().strip())

                if self.bond.text().strip() != '0':
                    bond = float(self.bond.text().strip())
                if self.angle.text().strip() != '0':
                    angle = float(self.angle.text().strip())
                if self.dihedral.text().strip() != '0':
                    dihedral = float(self.dihedral.text().strip())
                if self.improper.text().strip() != '0':
                    improper = float(self.improper.text().strip())
                if self.soft_sphere.text().strip() != '0':
                    soft_sphere = float(self.soft_sphere.text().strip())
                if self.lennard_jones.text().strip() != '0':
                    lennard_jones = float(self.lennard_jones.text().strip())
                if self.coulomb.text().strip() != '0':
                    coulomb = float(self.coulomb.text().strip())
                if self.h_bond.text().strip() != '0':
                    h_bond = float(self.h_bond.text().strip())
                if self.ca_distance.text().strip() != '0':
                    ca_distance = float(self.ca_distance.text().strip())
                if self.n_o_distance.text().strip() != '0':
                    n_o_distance = float(self.n_o_distance.text().strip())
                if self.phi_dihedral.text().strip() != '0':
                    phi_dihedral = float(self.phi_dihedral.text().strip())
                if self.psi_dihedral.text().strip() != '0':
                    psi_dihedral = float(self.psi_dihedral.text().strip())
                if self.omega_dihedral.text().strip() != '0':
                    omega_dihedral = float(self.omega_dihedral.text().strip())
                if self.chi1_dihedral.text().strip() != '0':
                    chi1_dihedral = float(self.chi1_dihedral.text().strip())
                if self.chi2_dihedral.text().strip() != '0':
                    chi2_dihedral = float(self.chi2_dihedral.text().strip())
                if self.chi3_dihedral.text().strip() != '0':
                    chi3_dihedral = float(self.chi3_dihedral.text().strip())
                if self.chi4_dihedral.text().strip() != '0':
                    chi4_dihedral = float(self.chi4_dihedral.text().strip())
                if self.disulfide_distance.text().strip() != '0':
                    disulfide_distance = float(self.disulfide_distance.text().strip())
                if self.disulfide_angle.text().strip() != '0':
                    disulfide_angle = float(self.disulfide_angle.text().strip())
                if self.disulfide_dihedral.text().strip() != '0':
                    disulfide_dihedral = float(self.disulfide_dihedral.text().strip())
                if self.lower_distance.text().strip() != '0':
                    lower_distance = float(self.lower_distance.text().strip())
                if self.upper_distance.text().strip() != '0':
                    upper_distance = float(self.upper_distance.text().strip())
                if self.sd_mn_distance.text().strip() != '0':
                    sd_mn_distance = float(self.sd_mn_distance.text().strip())
                if self.chi5_dihedral.text().strip() != '0':
                    chi5_dihedral = float(self.chi5_dihedral.text().strip())
                if self.phi_psi_dihedral.text().strip() != '0':
                    phi_psi_dihedral = float(self.phi_psi_dihedral.text().strip())
                if self.sd_sd_distance.text().strip() != '0':
                    sd_sd_distance = float(self.sd_sd_distance.text().strip())
                if self.xy_distance.text().strip() != '0':
                    xy_distance = float(self.xy_distance.text().strip())
                if self.nmr_distance.text().strip() != '0':
                    nmr_distance = float(self.nmr_distance.text().strip())
                if self.nmr_distance2.text().strip() != '0':
                    nmr_distance2 = float(self.nmr_distance2.text().strip())
                if self.min_distance.text().strip() != '0':
                    min_distance = float(self.min_distance.text().strip())
                if self.accessibility.text().strip() != '0':
                    accessibility = float(self.accessibility.text().strip())
                if self.density.text().strip() != '0':
                    density = float(self.density.text().strip())
                if self.absposition.text().strip() != '0':
                    absposition = float(self.absposition.text().strip())
                if self.dihedral_diff.text().strip() != '0':
                    dihedral_diff = float(self.dihedral_diff.text().strip())
                if self.SAXS.text().strip() != '0':
                    saxs = float(self.SAXS.text().strip())
                if self.symmetry.text().strip() != '0':
                    symmetry = float(self.symmetry.text().strip())

                env = environ()
                env.schedule_scale = physical.values(bond=bond, angle=angle, dihedral=dihedral, improper=improper,
                                                     soft_sphere=soft_sphere, lennard_jones=lennard_jones, coulomb=coulomb,
                                                     h_bond=h_bond, ca_distance=ca_distance, n_o_distance=n_o_distance,
                                                     phi_dihedral=phi_dihedral, psi_dihedral=psi_dihedral, omega_dihedral=omega_dihedral,
                                                     chi1_dihedral=chi1_dihedral, chi2_dihedral=chi2_dihedral, chi3_dihedral=chi3_dihedral,
                                                     chi4_dihedral=chi4_dihedral, disulfide_distance=disulfide_distance,
                                                     disulfide_angle=disulfide_angle, disulfide_dihedral=disulfide_dihedral,
                                                     lower_distance=lower_distance, upper_distance=upper_distance, sd_mn_distance=sd_mn_distance,
                                                     chi5_dihedral=chi5_dihedral, phi_psi_dihedral=phi_psi_dihedral, sd_sd_distance=sd_sd_distance,
                                                     xy_distance=xy_distance, nmr_distance=nmr_distance, nmr_distance2=nmr_distance2,
                                                     min_distance=min_distance, nonbond_spline=nonbond_spline, accessibility=accessibility,
                                                     density=density, absposition=absposition, dihedral_diff=dihedral_diff, saxs=saxs, symmetry=symmetry)
                env.io.atom_files_directory = ['../atom_files']
                env.edat.dynamic_sphere = True
                env.libs.topology.read(file='$(LIB)/top_heav.lib')
                env.libs.parameters.read(file='$(LIB)/par.lib')

                # get codes
                addL = []
                selected = self.model.selectedItems()
                for i in selected:
                    addL.append(i.text())
                for i in addL:
                    code = os.path.splitext(i)[0]
                    mdl = complete_pdb(env, os.path.join(self.path, i))
                    mdl.write(file=os.path.join(self.path, code + '.ini'))
                    # Select all atoms:
                    atmsel = selection(mdl)
                    # add restraint from file
                    if self.restraint.currentText() != 'Select':
                        mdl.restraints.append(file=os.path.join(self.path, self.restraint.currentText()))

                    mpdf = atmsel.energy()

                    # Create optimizer objects and set defaults for all further optimizations
                    cg = conjugate_gradients(output='REPORT', min_atom_shift=float(self.minAtomShift.text().strip()),
                                             residue_span_range=(int(self.From.text().strip()), int(self.To.text().strip())))
                    qn = quasi_newton(output='REPORT', min_atom_shift=float(self.minAtomShift.text().strip()),
                                      max_atom_shift=float(self.maxAtomShift.text().strip()),
                                      residue_span_range=(int(self.From.text().strip()), int(self.To.text().strip())))
                    md = molecular_dynamics(output='REPORT', cap_atom_shift=float(self.capAtomShift.text().strip()), md_time_step=float(self.mdTimeStep.text().strip()),
                                            init_velocities=True, md_return='FINAL', equilibrate=999999, guide_factor=0.0, guide_time=0.0, friction=0.0,
                                            residue_span_range=(int(self.From.text().strip()), int(self.To.text().strip())))

                    # Open a file to get basic stats on each optimization
                    trcfil = open(os.path.join(self.path, code + '_trace'), 'w')
                    # Run CG on the all-atom selection; write stats every 5 steps
                    if self.optimizeMethod.currentText() == 'Conjugate gradients' or self.optimizeMethod.currentText() == '1 & 3':
                        cg.optimize(atmsel,  max_iterations=int(self.maxIter.text().strip()), actions=actions.trace(5, trcfil))

                    if self.optimizeMethod.currentText() == 'Quasi newton':
                        qn.optimize(atmsel, max_iterations=int(self.maxIter.text().strip()), actions=actions.trace(5, trcfil))

                    if self.optimizeMethod.currentText() == 'Molecular dynamics' or self.optimizeMethod.currentText() == '1 & 3':
                        md.optimize(atmsel, temperature=float(self.temperature.text().strip()), max_iterations=int(self.maxIter.text().strip()),
                                    actions=[actions.trace(10, trcfil)])

                    mpdf = atmsel.energy()
                    mdl.write(file=os.path.join(self.path, code + '_optimized.pdb'))
                self.model.clear()
                self.model.addItems([i for i in os.listdir(self.path) if i.endswith('.pdb')])
                self.msg_1.setStyleSheet('color: green')
                self.msg_1.setText('Finished')
                self.optBackBoneBut.setEnabled(True)
        except ModuleNotFoundError:
            self.msg_1.setStyleSheet('color: red')
            self.msg_1.setText('MODELLER not found')
            self.optBackBoneBut.setEnabled(True)
        except Exception as er:
            self.msg_1.setStyleSheet('color: red')
            self.msg_1.setText('Error')
            self.optBackBoneBut.setEnabled(True)
            print(er)
            #print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))


# def main():
#     app = QApplication(sys.argv)
#     form = OptimizeModel()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()