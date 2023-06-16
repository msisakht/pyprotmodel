import os
import sys
import re
import json
import random
import threading
import time
import requests
import collections
import clipboard
import pandas as pd
import PyQt5
from PyQt5 import QtCore
from PyQt5.QtCore import QObject
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QMainWindow, QApplication, QHeaderView, QMenu
import tpl_search_seq
import config
import pd_model


class SeqSearch(QMainWindow, QObject, tpl_search_seq.Ui_Form):
    dir_signal = QtCore.pyqtSignal(str)
    # create info file
    os.makedirs('../templates', exist_ok=True)
    os.makedirs('../temp', exist_ok=True)
    if not os.path.exists(os.path.join(os.getcwd(), 'info')):
        with open('info', 'w') as f:
            json.dump(
                {'Path': os.path.join(os.getcwd(), '../templates'), 'Browsed_path': [],
                 'Target': '', 'multi_chain': '', 'add_ligands': '', 'add_ligand_rsrFile': '', 'add_ligand_feature_mat': '',
                 'use_symmetry': '', 'loop_res_num': '', 'symmetry': ''}, f)
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']
    # check for valid path
    infoFile['Browsed_path'] = [pth for pth in infoFile['Browsed_path'] if os.path.exists(pth)]

    search_result = {}
    selectedL = []
    colors_10 = ['#afd8aa', '#e88a5c', '#a8bdec', '#cca2c5', '#edd269', '#b4e4f7', '#ecbe82', '#ecc9c1', '#efbdc9', '#fbdce2',
                 '#d5e4b3', '#a5aaf7', '#8faaca', '#c7a7a7', '#b5c2d6', '#b4f7ec', '#89adb5', '#b497b8', '#bd9dd0', '#c2d697',
                 '#8dd6f5', '#bddbec', '#b6decc', '#a68cb2', '#94c191', '#b4c6f2', '#e7daa8', '#e0c834', '#d0d9df', '#dbc9e4']
    structure_dic = collections.OrderedDict({'Accession': [], 'Structure': [], 'Method': [], 'Resolution': [],
                                             'Description': []})
    model_dic = collections.OrderedDict({'Accession': [], 'Template': [], 'Method': [], 'Resolution': [], 'Similarity': [],
                                         'Identity': [], 'QMEAN': [], 'Q_Norm': [], 'Length': [], 'From': [], 'To': []})
    pdb_model_dic = collections.OrderedDict({'Accession': [], 'Method': [], 'Description': []})
    no_structureL = set()
    errorL = []

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        # working dir
        # add subdir of Templates dir
        self.dirL = [i for i in os.listdir('../templates') if os.path.isdir(os.path.join('../templates', i))]
        path = []
        for i in self.infoFile['Browsed_path']:
            if os.path.basename(i) != '':
                path.append(os.path.basename(i))
            else:
                path.append(i)
        self.workingDir.addItems(self.dirL + path)
        # set last dir to working dir
        if os.path.basename(self.path) != '':
            self.workingDir.setCurrentText(os.path.basename(self.path))
        else:
            self.workingDir.setCurrentText(self.path)
        self.workingDir.currentIndexChanged.connect(self.set_working_dir)
        self.workingDirShow.setText(self.path)
        self.browsDir.released.connect(self.browse_dir)
        # set modeller path
        self.connect_modeller()
        self.connectBut.released.connect(self.connect_modeller)
        # search type
        self.searchType.addItems(['Keywords', 'Accession numbers'])
        # return search num
        self.maxSearch.addItems([str(i) for i in range(5, 55, 5)])
        # search
        self.searchBut.released.connect(self.search_uniprot_th)
        self.search.returnPressed.connect(self.search_uniprot_th)
        # context menu
        self.accResult.installEventFilter(self)
        # display sequence
        self.seq_win.hide()
        # add items to check model
        self.accResult.clicked.connect(self.add_to_check_structure)
        # check models
        self.display_3d.addItems(['Experimental structure', 'PDB model', 'Model'])
        self.display_3d.currentIndexChanged.connect(self.check_display_3d)
        self.checkBut.released.connect(self.check_structure_th)
        self.set_working_dir()

    def set_working_dir(self):
        with open('info', 'r+') as f:
            jFile = json.load(f)
            path = jFile['Path']
            if self.workingDir.currentText() == 'templates':
                path = os.path.join(os.getcwd(), 'templates')
            if os.path.isdir(os.path.join('templates', self.workingDir.currentText())):
                path = os.path.join('templates', self.workingDir.currentText())
            else:
                for pth in jFile['Browsed_path']:
                    if self.workingDir.currentText() == os.path.basename(pth):
                        path = pth
            jFile['Path'] = os.path.abspath(path)
            f.seek(0)
            f.write(json.dumps(jFile))
            f.truncate()
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']
        self.workingDirShow.setText(path)
        self.dir_signal.emit(path)

    def connect_modeller(self):
        def read_path_from_reg():
            try:
                modeller_Path = os.path.abspath(os.path.join(config.Config.get_modeller_path().strip(), 'modlib'))
                self.modellerPath.setText(config.Config.get_modeller_path().strip())
                sys.path.insert(0, modeller_Path)
                from modeller import info
                config.Config.set_modeller_path(self.modellerPath.text().strip(), info.version)
                self.modellerConLable.setStyleSheet('color: green')
                self.modellerConLable.setText('Connected')
            except:
                self.modellerConLable.setStyleSheet('color: red')
                self.modellerConLable.setText('Not connected')
                pass

        if os.path.exists(self.modellerPath.text().strip()):
            modeller_Path = self.modellerPath.text().strip()
            sys.path.append(os.path.join(modeller_Path, 'modlib'))
            try:
                from modeller import info
                config.Config.set_modeller_path(modeller_Path, info.version)
                self.modellerConLable.setStyleSheet('color: green')
                self.modellerConLable.setText('Connected')
            except ModuleNotFoundError as er:
                print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
                read_path_from_reg()
        else:
            read_path_from_reg()

    def browse_dir(self):
        setPath = str(PyQt5.QtWidgets.QFileDialog.getExistingDirectory())
        with open('info', 'r+') as f:
            jFile = json.load(f)
            if not setPath:
                path = jFile['Path']
            else:
                path = os.path.abspath(setPath)
                if path not in jFile['Browsed_path']:
                    jFile['Browsed_path'].append(path)
            jFile['Path'] = os.path.abspath(path)
            f.seek(0)
            f.write(json.dumps(jFile))
            f.truncate()
        self.workingDir.clear()
        newPath = []
        for i in jFile['Browsed_path']:
            if os.path.basename(i) != '':
                newPath.append(os.path.basename(i))
            else:
                newPath.append(i)
        self.workingDir.addItems(['Templates'] + self.dirL + newPath)
        if os.path.basename(path) != '':
            self.workingDir.setCurrentText(os.path.basename(path))
        else:
            self.workingDir.setCurrentText(path)
        self.workingDirShow.setText(path)

    def eventFilter(self, source, event):
        try:
            if event.type() == QtCore.QEvent.ContextMenu and source is self.accResult:
                menu = QMenu()
                cp_acc = menu.addAction('Copy accession number(s)')
                show_seq = menu.addAction('Show sequence')
                action = menu.exec_(event.globalPos())
                if len(self.search_result) > 0:
                    if action == cp_acc:
                        self.copy_acc()
                        #item = source.itemAt(event.pos())
                        #print(item.text())
                else:
                    pass
                if action == show_seq:
                    self.seq_win.show()
                    self.show_seq()
        except:
            pass
        return super(SeqSearch, self).eventFilter(source, event)

    def copy_acc(self):
        addL = []
        copL = []
        selected = self.accResult.selectedItems()
        for i in selected:
            addL.append(i.text())
        for k, v in self.search_result.items():
            if k in addL:
                copL.append(v[0])
        clipboard.copy(' '.join(i for i in copL))

    def show_seq(self):
        addL = []
        selected = self.accResult.selectedItems()
        for i in selected:
            addL.append(i.text())
        for k, v in self.search_result.items():
            if k in addL:
                seq = ''.join(line.strip('\n') for line in v[2])
                self.seq_win.setWindowTitle('Sequence for: ' + v[0])
                self.seq_text.clear()
                self.seq_text.setText(seq)

    def search_uniprot_th(self):
        if len(self.search.text().strip()) != 0:
            threading.Thread(target=self.search_uniprot).start()
            self.searchBut.setEnabled(False)

    def search_uniprot(self):
        self.no_structureL.clear()
        self.structure_dic = collections.OrderedDict({'Accession': [], 'Structure': [], 'Method': [], 'Resolution': [],
                                                      'Description': []})
        self.model_dic = collections.OrderedDict({'Accession': [], 'Template': [], 'Method': [], 'Resolution': [],
                                                  'Similarity': [], 'Identity': [], 'QMEAN': [], 'Q_Norm': [],
                                                  'Length': [], 'From': [], 'To': []})
        try:
            self.searchBut.setText('Searchin...')
            s = re.findall(r'[A-Z\d\.]{4,}', self.search.text().upper())
            acL = list(set(s))
            if self.searchType.currentText() == 'Keywords':
                self.search_result.clear()
                self.selectedL.clear()
                self.run_search(self.search.text().strip())
                # check to get full data
                while True:
                    for v in self.search_result.values():
                        if len(v) < 3:
                            self.run_search(v[0])
                    break
                for k, v in self.search_result.items():
                    if len(v) == 2:
                        del self.search_result[k]

            elif self.searchType.currentText() == 'Accession numbers':
                self.search_result.clear()
                self.selectedL.clear()
                if len(acL) > 0:
                    for ac in acL:
                        self.run_search(ac)
                # check to get full data
                while True:
                    for v in self.search_result.values():
                        if len(v) < 3:
                            self.run_search(v[0])
                    break
                for k, v in self.search_result.items():
                    if len(v) == 2:
                        del self.search_result[k]
            if len(self.search_result) > 0:
                self.msgLabel.clear()
                self.accResult.clear()
                self.accResult.addItems([k for k in self.search_result.keys()])
                # add colors
                if len([k for k in self.search_result.keys()]) > len(acL) and any(i[4] != 1 for i in self.search_result.values()):
                    if len(acL) <= 10:
                        colors = self.colors_10
                    else:
                        colors = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(acL))]
                    for n, ac in enumerate(acL):
                        for i in [self.accResult.item(i) for i in range(self.accResult.count())]:
                            if i.text().startswith(ac):
                                i.setBackground(QColor(colors[n]))
            else:
                self.accResult.clear()
                self.msgLabel.setText('No result.')
            self.searchBut.setText('Search')
            self.searchBut.setEnabled(True)
        except Exception as er:
            print(er)
            print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
    
    def run_search(self, query):
        # connect to uniprot
        try:
            if len(self.search.text().strip()) >= 2:
                headers = {'Accept': 'application/json', }
                params = {'query': query, 'size': self.maxSearch.currentText().strip(),
                          'fields': 'accession,protein_name,length,organism_name,sequence', }
                res = requests.get('https://rest.uniprot.org/uniprotkb/search', headers=headers, params=params)
                resj = res.json()
                if len(resj['results']) > 0:
                    for n in range(len(resj['results'])):
                        accession = resj['results'][n]['primaryAccession']
                        name = resj['results'][n]['proteinDescription']['recommendedName']['fullName']['value']
                        length = resj['results'][n]['sequence']['length']
                        organism = resj['results'][n]['organism']['commonName']
                        sequence = resj['results'][n]['sequence']['value']
                        self.search_result[query + ': ' + accession + ' - ' + name.strip() + ' - (' +
                                    str(length) + ' AA) - (' + organism + ')'] = [accession,
                                    name.strip(), sequence.strip(), length, len(sequence)]
        except requests.exceptions.ConnectionError:
            self.accResult.clear()
            self.accResult.addItems(['Connection error.'])
        except Exception as er:
            self.accResult.clear()
            self.accResult.addItems([er])
            # print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            print(er)

    def add_to_check_structure(self):
        selected_idsL = [re.search(r'(\w+\d+)\s*-\s*.*', i.text()).group(1) for i in self.accResult.selectedItems()][:20]
        self.group_check_structure.setTitle(('Check structure For: ' + ', '.join(i for i in selected_idsL)))

    def check_structure_th(self):
        t = threading.Thread(target=self.check_structure)
        t.start()
        self.checkBut.setEnabled(False)
        self.msg_label.setText('')

    def check_structure(self):
        self.errorL.clear()
        os.makedirs('../temp', exist_ok=True)
        self.no_struLabel.setText('')
        selected_ids = [re.search(r'(\w+\d+)\s*-\s*.*', i.text()).group(1) for i in self.accResult.selectedItems()]
        for n, id in enumerate(selected_ids):
            try:
                if id in selected_ids:
                    self.checkBut.setText('Checking %d/%d...' % (n + 1, len(selected_ids)))
                    QApplication.processEvents()
                    if os.path.exists(os.path.join('../temp', id + '.json')):
                        pass
                    elif not os.path.exists(os.path.join('../temp', id + '.json')):
                        self.get_swiss_model(id)
                # show model features
                self.get_structure(id)
            except:
                self.errorL.append(id)
        if len(self.errorL) > 0:
            self.checkBut.setText('Check structure (%d Error)' % (len(self.errorL)))
            self.checkBut.setToolTip('Error checking: ' + ', '.join(i for i in self.errorL))
        else:
            self.checkBut.setText('Check structure')
        self.checkBut.setEnabled(True)

    def get_swiss_model(self, id):
        try:
            url = 'https://swissmodel.expasy.org/repository/uniprot/' + id + '.json'
            req = requests.get(url)
            req.raise_for_status()
            with open(os.path.join('../temp', id + '.json'), 'wb') as f:
                f.write(req.content)
            self.msg_label.setText('')
        except requests.exceptions.ConnectionError:
            self.errorL.append(id)
            self.msg_label.setStyleSheet('color: red')
            self.msg_label.setText('Connection error')
            QApplication.processEvents()
        except:
            pass

    def check_display_3d(self):
        if self.display_3d.currentText() == 'Experimental structure':
            dic = self.structure_dic
        elif self.display_3d.currentText() == 'PDB model':
            dic = self.pdb_model_dic
        else:
            dic = self.model_dic
        if len(self.structure_dic['Structure']) or len(self.model_dic['Template']) > 0:
            self.result_df = pd.DataFrame.from_dict(dic, orient='index').transpose()
            qt_model = pd_model.PdModel(self.result_df)
            self.table_results.setModel(qt_model)
            self.table_results.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

    def get_structure(self, id):
        pdb_idsL = []
        to_url = []
        interval = 3
        from_db = 'UniProtKB_AC-ID'
        to_db = 'PDB'
        # get pdb ids from uniprot
        try:
            def submit_id_mapping(from_db, to_db, id):
                API_URL = "https://rest.uniprot.org"
                session = requests.Session()
                request = session.post(f"{API_URL}/idmapping/run", data={"from": from_db, "to": to_db, "ids": id},)
                request.raise_for_status()
                return request.json()["jobId"]

            def check_id_mapping_results_ready(job_id):
                API_URL = "https://rest.uniprot.org"
                session = requests.Session()
                while True:
                    request = session.get(f"{API_URL}/idmapping/status/{job_id}")
                    request.raise_for_status()
                    j = request.json()
                    if "jobStatus" in j:
                        if j["jobStatus"] == "RUNNING":
                            print(f"Retrying in {interval}s")
                            time.sleep(interval)
                        else:
                            raise Exception(j["jobStatus"])
                    else:
                        return bool(j["results"] or j["failedIds"])

            def get_id_mapping_results_link(job_id):
                API_URL = "https://rest.uniprot.org"
                session = requests.Session()
                url = f"{API_URL}/idmapping/details/{job_id}"
                request = session.get(url)
                request.raise_for_status()
                return request.json()["redirectURL"]

            def get_id_mapping_results(url):
                session = requests.Session()
                request = session.get(url)
                request.raise_for_status()
                return request.text.splitlines()

            job_id = submit_id_mapping(from_db, to_db, id)
            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                results = get_id_mapping_results(link)
                data = json.loads(results[0])
                results_dic = data["results"]
                pdb_idsL = [result["to"] for result in results_dic]
            # get method, resolution and title of templates from rcsb
            if os.path.exists(os.path.join('../temp', 'pdb_data')):
                pass
            else:
                with open(os.path.join('../temp', 'pdb_data'), 'w') as f:
                    json.dump({}, f)
            with open(os.path.join('../temp', 'pdb_data'), 'r+') as f:
                jFile = json.load(f)
                for i in pdb_idsL:
                    if jFile.get(i, 0) != 0:
                        pass
                    else:
                        jFile[i] = []

                # check for pdb data
                for k, v in jFile.items():
                    if k in pdb_idsL and len(v) == 0:
                        to_url.append(k)

                # connect to rcsb
                if len(to_url) > 0:
                    url = "https://data.rcsb.org/graphql"
                    query = """{
                      entries(entry_ids: [""" + ', '.join(['"' + i + '"' for i in to_url]) + """]) {
                        rcsb_entry_info {
                          experimental_method
                          resolution_combined
                        }
                        struct {
                          title
                        }
                      }
                    }"""
                    req = requests.post(url, json={'query': query})
                    req_results = json.loads(req.content)
                    pdb = {'title': '', 'expmethod': '', 'resolution': ''}
                    for n in range(len(to_url)):
                        pdb['title'] = req_results['data']['entries'][n]['struct']['title']
                        pdb['expmethod'] = req_results['data']['entries'][n]['rcsb_entry_info']['experimental_method']
                        resolution = req_results['data']['entries'][n]['rcsb_entry_info']['resolution_combined']
                        if resolution:
                            pdb['resolution'] = \
                            req_results['data']['entries'][n]['rcsb_entry_info']['resolution_combined'][0]
                        else:
                            pdb['resolution'] = ''
                        try:
                            for k, v in jFile.items():
                                if k == to_url[n]:
                                    jFile[k] = [pdb.get('expmethod', 'None'), pdb.get('resolution', 'None'), pdb.get('title', 'None')]
                        except:
                            pass
                f.seek(0)
                f.write(json.dumps(jFile))
                f.truncate()
            # check for models?
            if len(pdb_idsL) > 0:
                pdb_data = open(os.path.join('../temp', 'pdb_data'), 'r')
                j = json.load(pdb_data)
                for k, v in j.items():
                    if not 'model' in v[0].lower():
                        if k in pdb_idsL:
                            self.structure_dic['Accession'].append(id)
                            self.structure_dic['Structure'].append(k)
                            self.structure_dic['Method'].append(v[0])
                            self.structure_dic['Resolution'].append(v[1])
                            self.structure_dic['Description'].append(v[2])
                    else:
                        self.pdb_model_dic['Accession'].append('UniProt: ' + id + ', PDB: ' + k)
                        self.pdb_model_dic['Method'].append(v[0])
                        self.pdb_model_dic['Description'].append(v[2])

                if self.display_3d.currentText() == 'Experimental structure':
                    dic = self.structure_dic
                elif self.display_3d.currentText() == 'PDB model':
                    dic = self.pdb_model_dic
                else:
                    dic = self.model_dic
                # display results
                self.result_df = pd.DataFrame.from_dict(dic, orient='index').transpose()
                qt_model = pd_model.PdModel(self.result_df)
                self.table_results.setModel(qt_model)
                self.table_results.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            else:
                self.get_all_model_data(id)
        except requests.exceptions.ConnectionError:
            self.msg_label.setStyleSheet('color: red')
            self.msg_label.setText('Connection error')
            QApplication.processEvents()
        except:
            self.errorL.append(id)
            pass

    def get_all_model_data(self, id):
        for k, v in self.search_result.items():
            try:
                if v[0] == id and os.path.exists(os.path.join('../temp', v[0] + '.json')):
                    # ___ swiss model ___
                    identity = []
                    similarity = []
                    qmean = []
                    qmean_norm = []
                    template = []
                    frm = []
                    to = []
                    l = []
                    to_url = []
                    with open(os.path.join('../temp', v[0] + '.json'), 'r') as f:
                        j = json.load(f)
                        for n in range(len(j['result']['structures'])):
                            try:
                                identity.append(round(float(j['result']['structures'][n]['identity']), 2))
                            except:
                                identity.append('-')
                            try:
                                similarity.append(round(float(j['result']['structures'][n]['similarity']), 2))
                            except:
                                similarity.append('-')
                            try:
                                qmean.append(round(float(j['result']['structures'][n]['qmean']), 2))
                            except:
                                qmean.append('-')
                            try:
                                qmean_norm.append(round(float(j['result']['structures'][n]['qmean_norm']), 2))
                            except:
                                qmean_norm.append('-')
                            try:
                                template.append(j['result']['structures'][n]['template'])
                            except:
                                template.append('-')
                            try:
                                frm.append(j['result']['structures'][n]['from'])
                            except:
                                frm.append('-')
                            try:
                                to.append(j['result']['structures'][n]['to'])
                            except:
                                to.append('-')

                    # add pdb data to file
                    for i in template:
                        l.append(re.search(r'(\w{,5})', i).group(1))
                    if os.path.exists(os.path.join('../temp', 'pdb_data')):
                        pass
                    else:
                        with open(os.path.join('../temp', 'pdb_data'), 'w') as f:
                            json.dump({}, f)
                    with open(os.path.join('../temp', 'pdb_data'), 'r+') as f:
                        jFile = json.load(f)
                        for i in l:
                            if jFile.get(i, 0) != 0:
                                pass
                            else:
                                jFile[i] = []
                                to_url.append(i)
                        # connect to rcsb
                        if len(to_url) > 0:
                            url = "https://data.rcsb.org/graphql"
                            query = """{entries(entry_ids: [""" + ', '.join(['"' + i + '"' for i in to_url]) + """]) {
                                                    rcsb_entry_info {
                                                      experimental_method
                                                      resolution_combined
                                                    }
                                                    struct {
                                                      title
                                                    }
                                                  }
                                                }"""
                            req = requests.post(url, json={'query': query})
                            req_results = json.loads(req.content)
                            pdb = {'title': '', 'expmethod': '', 'resolution': ''}
                            for n in range(len(to_url)):
                                pdb['title'] = req_results['data']['entries'][n]['struct']['title']
                                pdb['expmethod'] = req_results['data']['entries'][n]['rcsb_entry_info']['experimental_method']
                                pdb['resolution'] = '-'
                                try:
                                    for k, v in jFile.items():
                                        if k == to_url[n]:
                                            jFile[k] = [pdb.get('expmethod', 'None'), pdb.get('resolution', 'None'), pdb.get('title', 'None')]
                                except:
                                    pass
                        f.seek(0)
                        f.write(json.dumps(jFile))
                        f.truncate()

                    # print data
                    if len(template) > 0:
                        pdb_data = open(os.path.join('../temp', 'pdb_data'), 'r')
                        j = json.load(pdb_data)
                        for n in range(len(template)):
                            self.model_dic['Accession'].append(id)
                            self.model_dic['Template'].append(template[n])
                            self.model_dic['Method'].append([value[0] for key, value in j.items() if key in template[n]][0])
                            self.model_dic['Resolution'].append([value[1] for key, value in j.items() if key in template[n]][0])
                            self.model_dic['Similarity'].append(similarity[n])
                            self.model_dic['Identity'].append(identity[n])
                            self.model_dic['QMEAN'].append(qmean[n])
                            self.model_dic['Q_Norm'].append(qmean_norm[n])
                            self.model_dic['Length'].append(v[3])
                            self.model_dic['From'].append(frm[n])
                            self.model_dic['To'].append(to[n])

                        pdb_data.close()
                        # show models
                        if self.display_3d.currentText() == 'Model':
                            dic = self.model_dic
                        elif self.display_3d.currentText() == 'PDB model':
                            dic = self.pdb_model_dic
                        else:
                            dic = self.structure_dic
                        # display results
                        self.result_df = pd.DataFrame.from_dict(dic, orient='index').transpose()
                        qt_model = pd_model.PdModel(self.result_df)
                        self.table_results.setModel(qt_model)
                        self.table_results.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
                    else:
                        self.no_structureL.add(id)
                        self.no_struLabel.setStyleSheet('color: red')
                        self.no_struLabel.setText('No structure for: ' + ', '.join(i for i in self.no_structureL))
            except IndexError:
                continue
            except Exception as er:
                print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
                pass


# def main():
#     app = QApplication(sys.argv)
#     form = SeqSearch()
#     form.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()