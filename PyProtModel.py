import os
import sys
import json
from winreg import *
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication
import tpl_tabs
import seq_search
import temp_search
import select_template
import repair_temp
import align_temp
import def_restraint
import add_ligand
import build_model
import optimize_model
import evaluate_model
import refine_loop


class PyProtModel(QMainWindow, tpl_tabs.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.setFixedSize(979, 605)
        #
        try:
            modelVer = self.get_modeller_version()
            self.setWindowTitle('PyProtModel V.1 (' + modelVer + ')')
        except:
            self.setWindowTitle('PyProtModel V.1')
            pass
        self.setFixedSize(self.size())
        #
        self.search_seq = QtWidgets.QWidget()
        self.search_seq.setObjectName("search_seq")
        self.tabWidget.addTab(seq_search.SeqSearch(), "Search Sequence")

        self.search_temp = QtWidgets.QWidget()
        self.search_temp.setObjectName("search_temp")
        self.tabWidget.addTab(temp_search.TempSearch(), "Search Templates")

        self.select_temp = QtWidgets.QWidget()
        self.select_temp.setObjectName("select_temp")
        self.tabWidget.addTab(select_template.SelectTemp(), "Select Templates")

        self.repair_missing_atom = QtWidgets.QWidget()
        self.repair_missing_atom.setObjectName("repair_missing_atom")
        self.tabWidget.addTab(repair_temp.RepairTemp(), "Repair Templates")

        self.align_templates = QtWidgets.QWidget()
        self.align_templates.setObjectName("align_templates")
        self.tabWidget.addTab(align_temp.AlignTemp(), "Align Templates")

        self.def_restraint = QtWidgets.QWidget()
        self.def_restraint.setObjectName("def_restraint")
        self.tabWidget.addTab(def_restraint.DefineRestraint(), "Define Restraint")

        self.add_lig = QtWidgets.QWidget()
        self.add_lig.setObjectName("add_lig")
        self.tabWidget.addTab(add_ligand.AddLigand(), "Add Ligand")

        self.build_model = QtWidgets.QWidget()
        self.build_model.setObjectName("build_model")
        self.tabWidget.addTab(build_model.BuildModel(), "Build Model")

        self.opt_model = QtWidgets.QWidget()
        self.opt_model.setObjectName("opt_model")
        self.tabWidget.addTab(optimize_model.OptimizeModel(), "Optimize Model")

        self.loop_ref = QtWidgets.QWidget()
        self.loop_ref.setObjectName("loop_ref")
        self.tabWidget.addTab(refine_loop.RefineLoop(), "Refine Model")

        self.eval_model = QtWidgets.QWidget()
        self.eval_model.setObjectName("eval_model")
        self.tabWidget.addTab(evaluate_model.EvaluateModel(), "Evaluate Model")

    def get_modeller_version(self):
        try:
            pyprotmodel_key = OpenKey(HKEY_CURRENT_USER, r'SOFTWARE\PyProtModel', 0, KEY_READ)
            [pathVal, regtype] = (QueryValueEx(pyprotmodel_key, 'MODELLER_VERSION'))
            CloseKey(pyprotmodel_key)
            return pathVal
        except:
            pass


def welcome():
    os.system('color 4f')
    os.system('cls')
    print('\nDEPARTMENT OF MEDICAL BIOCHEMISTRY, FACULTY OF MEDICINE, SHIRAZ, IRAN\n\n')
    print('PyProtModel Version 1.0')
    print(' _____   __     __  __  __    ____    _____    ______   _      ')
    print('|  __ \  \ \   / / |  \/  |  / __ \  |  __ \  |  ____| | |     ')
    print('| |__) |  \ \_/ /  | \  / | | |  | | | |  | | | |__    | |     ')
    print('|  ___/    \   /   | |\/| | | |  | | | |  | | |  __|   | |     ')
    print('| |         | |    | |  | | | |__| | | |__| | | |____  | |____ ')
    print('|_|         |_|    |_|  |_|  \____/  |_____/  |______| |______|\n')
    print('\nFor any help and support please send an email to: mohsen.sisakht@gmail.com\n')
    os.system('pause')
    os.system('cls')


def main():
    global app
    welcome()
    os.system('color 0f')
    app = QApplication(sys.argv)
    form = PyProtModel()
    form.show()
    app.exec_()
    os.system('pause')


if __name__ == '__main__':
    main()
