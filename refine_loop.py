import os
import sys
import json
import shutil
import threading
import pandas as pd
from Bio.PDB import PDBParser
from PyQt5.QtWidgets import QMainWindow, QApplication, QHeaderView, QFileDialog
import pd_model
import tpl_refine_loop
import config


class RefineLoop(QMainWindow, tpl_refine_loop.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    pdb_fileL = []
    rsr_fileL = []
    browse_pdb = False
    refLevelL = ['Very fast', 'Fast', 'Slow', 'Very slow', 'Slow large']
    addedResChainDic = {}
    assessL = ['GA341', 'DOPE', 'DOPE-HR', 'Normalized DOPE', 'SOAP_PP', 'SOAP-Loop', 'SOAP-Peptide', 'SOAP-Protein']

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.model.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.pdb')])
        self.model.currentIndexChanged.connect(self.get_chain)
        self.browseModel.released.connect(self.browse_model_th)
        #
        self.chainFrom.currentIndexChanged.connect(self.get_res_num)
        self.chainTo.currentIndexChanged.connect(self.get_res_num)
        #
        self.addButton.released.connect(self.add_res_chain)
        self.minButton.released.connect(self.min_res_chain)
        #
        self.refLevel.addItems(self.refLevelL)
        #
        self.numModel.setText('50')
        #
        self.nameModel.setText('refined')
        #
        self.restraint.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.rsr')])
        #
        self.assess.addItems(self.assessL)
        #
        self.group_loop_result.hide()
        #
        self.pdb_fileL = [i for i in os.listdir(self.path) if i.endswith('.pdb')]
        self.rsr_fileL = [i for i in os.listdir(self.path) if i.endswith('.rsr')]
        #
        self.saveCsvBut.released.connect(self.save_result_csv)
        #
        self.refineBut.released.connect(self.refine_loop_th)

    def enterEvent(self, QEvent):
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']
        self.path = path
        #
        pdb_files = []
        rsr_files = []
        for i in os.listdir(self.path):
            if i.endswith('.pdb') and i not in pdb_files:
                pdb_files.append(i)
            if i.endswith('.rsr') and i not in rsr_files:
                rsr_files.append(i)
        while self.pdb_fileL != pdb_files:
            if not self.browse_pdb:
                self.model.clear()
                self.model.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.pdb')])
                self.pdb_fileL = pdb_files
        while self.rsr_fileL != rsr_files:
            self.restraint.clear()
            self.restraint.addItems(['Select'] + [i for i in os.listdir(path) if i.endswith('.rsr')])
            self.rsr_fileL = rsr_files

    def browse_model_th(self):
        threading.Thread(target=self.browse_model).start()

    def browse_model(self):
        user_pdb = QFileDialog.getOpenFileNames()[0]
        if user_pdb:
            for file in user_pdb:
                shutil.copy(file, self.path)
            self.model.clear()
            self.model.addItems([os.path.basename(i) for i in os.listdir(self.path) if i.endswith('.pdb')])
            self.model.setCurrentText(os.path.basename(user_pdb[0]))
            self.browse_pdb = True

    def get_chain(self):
        chainL = []
        numFL = []
        numTL = []
        self.chainFrom.clear()
        self.chainTo.clear()
        self.resFrom.clear()
        self.resTo.clear()
        try:
            if self.model.currentText() != 'Select':
                p = PDBParser(PERMISSIVE=1)
                s = p.get_structure(self.model.currentText(), os.path.join(self.path, self.model.currentText()))
                # get chains
                for chain in s.get_chains():
                    chainL.append(chain.id)
                if len(chainL) > 0:
                    self.chainFrom.addItems(chainL)
                    self.chainTo.addItems(chainL)
                else:
                    self.chainFrom.addItems([''])
                    self.chainTo.addItems([''])
                # get res number of each chain
                for num in s[0][self.chainFrom.currentText()]:
                    numFL.append(num.id[1])
                for num in s[0][self.chainTo.currentText()]:
                    numTL.append(num.id[1])
                if len(numFL) > 0:
                    self.resFrom.addItems(str(i) for i in numFL)
                    self.resTo.addItems(str(i) for i in numTL)
                else:
                    self.resFrom.addItems([''])
                    self.resTo.addItems([''])
            else:
                self.chainFrom.addItems([''])
                self.chainTo.addItems([''])
                self.resFrom.addItems([''])
                self.resTo.addItems([''])
        except Exception as er:
            # print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            pass

    def get_res_num(self):
        numFL = []
        numTL = []
        self.resFrom.clear()
        self.resTo.clear()
        #
        if self.chainFrom.currentText() != '':
            try:
                p = PDBParser(PERMISSIVE=1)
                s = p.get_structure(self.model.currentText(), os.path.join(self.path, self.model.currentText()))
                # get res number of each chain
                for num in s[0][self.chainFrom.currentText()]:
                    numFL.append(num.id[1])
                for num in s[0][self.chainTo.currentText()]:
                    numTL.append(num.id[1])
                if len(numFL) > 0:
                    self.resFrom.addItems(str(i) for i in numFL)
                else:
                    self.resFrom.addItems([''])
                if len(numTL) > 0:
                    self.resTo.addItems(str(i) for i in numTL)
                else:
                    self.resTo.addItems([''])
            except Exception as er:
                # print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
                pass

    def add_res_chain(self):
        k = [self.model.currentText(), self.chainFrom.currentText(), self.resFrom.currentText(), self.chainTo.currentText(), self.resTo.currentText()]
        # add data to dic
        if len(k[2].strip()) > 0 and len(k[4].strip()) > 0:
            self.addedResChainDic[k[0] + ', ' + k[1] + ': ' + k[2] + ', ' + k[3] + ': ' + k[4]] = [k[0], k[1], k[2], k[3], k[4]]
        dic = dict(self.addedResChainDic)
        for key, value in dic.items():
            if value[0] != k[0]:
                del self.addedResChainDic[key]
        self.resNum.clear()
        self.resNum.addItems([k for k in self.addedResChainDic.keys()])
        # write to json info file
        with open('info', 'r+') as f:
            jFile = json.load(f)
            jFile['loop_res_num'] = [v for v in self.addedResChainDic.values()]
            f.seek(0)
            f.write(json.dumps(jFile))
            f.truncate()

    def min_res_chain(self):
        selected = self.resNum.selectedItems()
        for i in selected:
            del self.addedResChainDic[i.text()]
        self.resNum.clear()
        self.resNum.addItems([k for k in self.addedResChainDic.keys()])
        # update json info file
        with open('info', 'r+') as f:
            jFile = json.load(f)
            jFile['loop_res_num'] = [v for v in self.addedResChainDic.values()]
            f.seek(0)
            f.write(json.dumps(jFile))
            f.truncate()

    def refine_loop_th(self):
        threading.Thread(target=self.refine_loop).start()
        self.refineBut.setEnabled(False)

    def refine_loop(self):
        self.msg.setText('Processing...')
        modeller_path = config.Config.get_modeller_path()
        sys.path.insert(0, modeller_path)
        fileL_before = [i for i in os.listdir(os.getcwd())]
        try:
            from modeller import log, environ, soap_protein_od, soap_pp, soap_loop, soap_peptide, restraints
            from modeller.automodel import assess, refine, loopmodel
            assessFuncL = (
            assess.GA341, assess.DOPE, assess.DOPEHR, assess.normalized_dope, soap_pp, soap_loop, soap_peptide,
            soap_protein_od)
            log.verbose()
            env = environ()
            env.io.atom_files_directory = './:../atom_files'
            env.libs.topology.read(file='$(LIB)/top_heav.lib')
            env.libs.parameters.read(file='$(LIB)/par.lib')
            # get refinemebt level
            md_level = refine.very_fast
            if self.refLevel.currentText() == 'Very fast':
                md_level = refine.very_fast
            elif self.refLevel.currentText() == 'Fast':
                md_level = refine.fast
            elif self.refLevel.currentText() == 'Slow':
                md_level = refine.slow
            elif self.refLevel.currentText() == 'Very slow':
                md_level = refine.very_slow
            elif self.refLevel.currentText() == 'Slow large':
                md_level = refine.slow_large

            # get csr file
            if self.restraint.currentText() == 'Select':
                restraint = None
            else:
                restraint = open(os.path.join(self.path, self.restraint.currentText()), 'r')

            # gat assess method
            getAssessL = []
            selection = self.assess.selectedItems()
            for i in selection:
                getAssessL.append(i.text())
            idx = []
            for i in self.assessL:
                if i in getAssessL:
                    idx.append(self.assessL.index(i))
            assess_methods = ([i for i in assessFuncL if assessFuncL.index(i) in idx])
            # run refinement

            class RefineLoop(loopmodel):

                def select_loop_atoms(self):
                    file = open('info')
                    infoFile = json.load(file)
                    sys.path.insert(0, modeller_path)
                    from modeller import selection

                    res_ranges = infoFile['loop_res_num']
                    return selection(
                        self.residue_range(num, chain) for num, chain in [(res_ranges[n][2] + ':' + res_ranges[n][1],
                                                                           res_ranges[n][4] + ':' + res_ranges[n][3])
                                                                          for n in range(len(res_ranges))])
                    #return selection(self.atom_range('CA:1:A', 'CB:10:A'))
                    #return selection(self.residue_range('1:A', '10:A'))

                def special_restraints(self, aln):
                    rsr = self.restraints
                    try:
                        rsr.append(restraint)
                    except:
                        pass

            r_loop = RefineLoop(env, inimodel=os.path.join(self.path, self.model.currentText()), sequence=self.nameModel.text().strip(),
                                loop_assess_methods=assess_methods)

            r_loop.loop.starting_model = 1
            r_loop.loop.ending_model = int(self.numModel.text().strip())
            r_loop.loop.md_level = md_level
            r_loop.make()
            fileL_after = [i for i in os.listdir(os.getcwd())]
            # get model files
            for i in fileL_after:
                if i not in fileL_before:
                    if os.path.exists(os.path.join(self.path, i)):
                        os.unlink(os.path.join(self.path, i))
                    shutil.move(os.path.join(os.getcwd(), i), self.path)
            self.msg.setStyleSheet('color: green')
            self.msg.setText('Finished')
            self.refineBut.setEnabled(True)
            # display results
            self.group_loop_result.show()
            # set tableView in full width
            self.table_loop_result.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            #
            model_output = r_loop.loop.outputs
            for n in range(len(r_loop.loop.outputs)):
                if 'GA341 score' in r_loop.loop.outputs[n]:
                    model_output[n]['GA341 score'] = r_loop.loop.outputs[n]['GA341 score'][0]
            self.model_df = pd.DataFrame(model_output)
            self.model_df.drop(['failure', 'num', 'loopnum', 'pdfterms'], axis=1, inplace=True)
            self.model_df = self.model_df[self.model_df.columns[::-1]]
            qt_model = pd_model.PdModel(self.model_df)
            self.table_loop_result.setModel(qt_model)
        except ModuleNotFoundError:
            self.msg.setStyleSheet('color: red')
            self.msg.setText('MODELLER not found')
            self.refineBut.setEnabled(True)
        except Exception as er:
            self.msg.setStyleSheet('color: red')
            self.msg.setText('Error')
            self.refineBut.setEnabled(True)
            # print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            print(er)

    def save_result_csv(self):
        try:
            self.msg_csv.setStyleSheet('color: #000000')
            self.msg_csv.setText('Saved')
            QApplication.processEvents()
            pd.DataFrame.to_csv(self.model_df, os.path.join(self.path, 'model_refinement.csv'), index=False)
        except Exception as er:
            self.msg_csv.setStyleSheet('color: red')
            self.msg_csv.setText('Error')
            print(er)
            pass


# def main():
#     app = QApplication(sys.argv)
#     form = RefineLoop()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()