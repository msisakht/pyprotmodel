import os
import sys
import re
import json
import csv
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import threading
import requests
import fileinput
import bs4
from Bio import SearchIO
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBList
from PyQt5.QtWidgets import QMainWindow, QApplication, QCheckBox, QHBoxLayout, QVBoxLayout, QHeaderView
import tpl_select_temp
import pd_model
import run_clustalo


class SelectTemp(QMainWindow, tpl_select_temp.Ui_Form):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']

    blast_fileL = []
    pdbL = []
    c_id = {}
    pdb_index = set()

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        # check for blast files
        try:
            self.blast_fileL = [i for i in os.listdir(self.path) if i.endswith('.blast')]
            self.blastFile.addItems([i for i in os.listdir(self.path) if
                                     os.path.isfile(os.path.join(self.path, i)) and i.endswith('.blast')])
        except:
            pass
        # blast files
        self.blastFile.activated.connect(self.get_templates)
        # custom ids
        self.pdb_selected.textEdited.connect(self.custom_id)
        # show templates
        self.showBut.released.connect(self.get_templates)
        self.FBRec.returnPressed.connect(self.get_templates)
        # download templates
        self.downloadBut.released.connect(self.download_pdb_th)
        # compare templates
        self.compareBut.released.connect(self.compare_template_th)
        # view compare result
        self.viewBut.setVisible(False)
        self.viewBut.released.connect(self.view_compare_result)

    def enterEvent(self, QEvent):
        file = open('info')
        infoFile = json.load(file)
        path = infoFile['Path']

        while self.path != path:
            self.blastFile.clear()
            self.blastFile.addItems([i for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and i.endswith('.blast')])
            self.path = path
        #
        blast_files = []
        for i in os.listdir(self.path):
            if i.endswith('.blast') and i not in blast_files:
                blast_files.append(i)
        while self.blast_fileL != blast_files:
            self.blastFile.clear()
            self.blastFile.addItems([i for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and i.endswith('.blast')])
            self.blast_fileL = blast_files

    def get_templates(self):
        try:
            # get accession num of primary protein
            for i in os.listdir(self.path):
                if i.endswith('_info.txt'):
                    with open(os.path.join(self.path, i), 'r') as f:
                        self.s = re.search(r'Accession: (.*)\nProgram: (.*)', f.read())
            self.programTypeLabel.setText(self.s.group(2))
            # ___searchIO read file___
            if self.s.group(2) == 'NCBI BLASTp' or self.s.group(2) == 'PSI-BLAST':
                blast_type = 'blast-text'
                read = SearchIO.read(os.path.join(self.path, self.blastFile.currentText()), 'blast-text')
                readNum = read[:int(self.FBRec.text())]
                self.show_result(blast_type, self.s.group(1), readNum)
            elif self.s.group(2) == 'FASTA BLAST Suite' or self.s.group(2) == 'PSI-Search':
                # ___bs4 read file___
                blast_type = 'xml'
                file = open(os.path.join(self.path, self.blastFile.currentText()))
                soup = bs4.BeautifulSoup(file, 'lxml')
                self.hit = soup.find_all('hit')
                readNum = self.hit[:int(self.FBRec.text())]
                self.show_result(blast_type, self.s.group(1), readNum)
        except requests.exceptions.ConnectionError:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Connection error')
            print('Connection error')
        except Exception as er:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Error')
            print(er)

    def show_result(self, blast_type, id, hit):
        # get pdb data
        l = []
        to_url = []
        for n in range(len(hit)):
            if blast_type == 'blast-text':
                l.append(re.search(r'(\w*:)?(\w{4})(_\w+)?', hit.hsps[n].hit_id[:18].strip()).group(2))
            elif blast_type == 'xml':
                l.append(re.search(r'(\w*:)?(\w{4})(_\w+)?', hit[n]['id'][:18].strip()).group(2))

        # get method and resolution of templates from rcsb
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
            # check for pdb data
            for k, v in jFile.items():
                if k in l[:int(self.FBRec.text())] and len(v) == 0:
                    to_url.append(k)
            # conncet to rcsb
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
                    resolution = req_results['data']['entries'][n]['rcsb_entry_info']['resolution_combined']
                    if resolution:
                        pdb['resolution'] = req_results['data']['entries'][n]['rcsb_entry_info']['resolution_combined'][0]
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

        header = ['Download', 'Template', 'Method', 'Resolution', 'Length/Range', 'Identity ', 'E value', 'Bit score',
                        'Description']
        templateData = []
        #
        pdb_data = open(os.path.join('../temp', 'pdb_data'), 'r')
        j = json.load(pdb_data)
        #
        if blast_type == 'blast-text':
            for n in range(len(hit.hsps)):
                try:
                    targLeng = hit.hits[n].seq_len
                    templateData.append(
                        [QCheckBox(''), hit.hsps[n].hit_id[:18], [v[0] for k, v in j.items() if k in hit.hsps[n].hit_id[:18]][0],
                         [v[1] for k, v in j.items() if k in hit.hsps[n].hit_id[:18]][0],
                         str(targLeng).strip() + ' / ' + str(hit.hits[n].fragments[0]._hit_range_get()),
                         hit.hsps[n].ident_num, hit.hsps[n].evalue, hit.hsps[n].bitscore,
                         re.search(r'mol:protein\s*length:\d*\s*(.*)', hit.hsps[n].hit_description).group(1)])
                    self.pdbL.append([n, re.search(r'(\w*:)?(\w{4})(_\w+)?', hit.hsps[n].hit_id[:18]).group(2),
                                      [v[1] for k, v in j.items() if k in hit.hsps[n].hit_id[:18]][0]])
                except IndexError:
                    targLeng = '-'
                    continue
        #
        elif blast_type == 'xml':
            # check id attr existanse
            try:
                for n in range(len(hit)):
                    if hit[n]['id']:
                        pass
            except:
                for n in range(len(hit)):
                    hit[n]['id'] = hit[n]['ac']
            #
            for n in range(len(hit)):
                templateData.append([QCheckBox(''), hit[n]['id'][:18], [v[0] for k, v in j.items() if k in hit[n]['id'][:18]][0],
                                    [v[1] for k, v in j.items() if k in hit[n]['id'][:18]][0], hit[n]['length'] + '/' +
                                    '(' + hit[n].alignment.matchseq['start'] + ', ' + hit[n].alignment.matchseq[
                                    'end'] + ')', hit[n].identity.getText(), hit[n].expectation.getText(),
                                    hit[n].bits.getText(), re.search(r'mol:protein\s*length:\d*\s*(.*)', hit[n]['description']).group(1)])
                self.pdbL.append([n, re.search(r'(\w*:)?(\w{4})(_\w+)?', hit[n]['id'][:18]).group(2),
                                  [v[1] for k, v in j.items() if k in hit[n]['id'][:18]][0]])

        pdb_data.close()
        # write csv file
        if self.saveCSV.isChecked():
            csvFile = open(os.path.join(self.path, id + '_templates.csv'), 'w', newline='')
            csvWriter = csv.writer(csvFile, delimiter=',', )
            csvWriter.writerow(
                ['Template', 'Method', 'Resolution', 'Length/Range', 'Identity ', 'E value', 'Bit score', 'Description'])
            for row in templateData:
                csvWriter.writerow(row[1:])
            csvFile.close()
        # show templates
        if len(templateData) > 0:
            model = pd_model.TableModel(templateData, header)
            self.template_table.setModel(model)
            layout = QHBoxLayout()
            Vlayout = QVBoxLayout()
            Vlayout.addWidget(self.template_table)
            Vlayout.addLayout(layout)
            self.setLayout(Vlayout)
            self.pdb_index.clear()
            # self.pdb_selected.clear()
            self.template_table.clicked.connect(self.get_index)
            self.template_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.msgLabel.setStyleSheet('color: #000000')
            self.msgLabel.setText('')
            QApplication.processEvents()
        else:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('No templates')
            QApplication.processEvents()

    def get_index(self, clicked):
        try:
            idx = clicked.row()
            data = clicked.data()
            if data.isChecked():
                self.pdb_index.add(clicked.row())
            else:
                if idx in self.pdb_index:
                    self.pdb_index.remove(idx)
            # display selected pdb cods
            all_selected = list(set(i[1] for i in self.pdbL if i[0] in self.pdb_index))
            pdb_selected_str = ', '.join(list(set(i for i in all_selected)) + self.custom_id())
            self.pdb_selected.setText(pdb_selected_str)
        except Exception as er:
            print(er)
            pass

    def custom_id(self):
        custom_ids = set()
        ids = [i for i in re.findall(r'(\w{4})', self.pdb_selected.text())]
        for i in ids:
            if i not in set(i[1] for i in self.pdbL):
                custom_ids.add(i)
        return list(custom_ids)

    def download_pdb_th(self):
        threading.Thread(target=self.download_pdb).start()

    def download_pdb(self):
        pdbdown = PDBList()
        # add entered pdb code
        u_code = re.findall(r'(\w{4})', self.pdb_selected.text())
        errorL = []
        if len(u_code) > 0:
            self.downloadBut.setEnabled(False)
            for n, i in enumerate(u_code):
                try:
                    if os.path.exists(os.path.join(self.path, i + '.pdb')):
                        pass
                    else:
                        self.msgLabel.setStyleSheet('color: black')
                        self.msgLabel.setText('Downloading %s (%d/%d)' % (i + '.pdb ', n + 1, len(u_code)))
                        pdbdown.retrieve_pdb_file(pdb_code=i, pdir=self.path, file_format='pdb', overwrite=True, obsolete=False)
                        os.rename(os.path.join(self.path, 'pdb' + i.lower() + '.ent'), os.path.join(self.path, i + '.pdb'))
                        if not os.path.exists(os.path.join(self.path, i + '.pdb')):
                            raise Exception
                except Exception as er:
                    errorL.append(i)
                    continue
        if len(errorL) == 0:
            self.msgLabel.setText('')
        else:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Error in downloading: ' + ', '.join(i for i in errorL))
        # get c_id resolution data
        if len(self.custom_id()) > 0:
            url = "https://data.rcsb.org/graphql"
            query = """{entries(entry_ids: [""" + ', '.join(['"' + i + '"' for i in self.custom_id()]) + """]) {
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
            for n, i in enumerate(self.custom_id()):
                resolution = req_results['data']['entries'][n]['rcsb_entry_info']['resolution_combined']
                if resolution:
                    pass
                else:
                    resolution = ['']
                self.pdbL.append([len(self.custom_id()) + n, i, resolution[0]])
        self.downloadBut.setEnabled(True)

    def compare_template_th(self):
        u_code = re.findall(r'(\w{4})', self.pdb_selected.text())
        if len(u_code) > 0:
            self.compareBut.setEnabled(False)
            self.compare_th = threading.Thread(target=self.compare_template)
            self.compare_th.start()

    def compare_template(self):
        self.download_pdb()
        try:
            self.msgLabel.setStyleSheet('color: black')
            self.msgLabel.setText('Comparing templates...')
            QApplication.processEvents()
            # get pdb file & chains
            p = PDBParser()
            for file in re.findall(r'(\w{4})', self.pdb_selected.text()):
                if os.path.exists(os.path.join(self.path, file + '.pdb')):
                    for i in self.pdbL:
                        if file.lower() == i[1].lower():
                            self.c_id[file] = str(i[2])
            # pdb to fasta
            seq = []
            for k in self.c_id.keys():
                s = p.get_structure(os.path.join(self.path, k + '.pdb'), os.path.join(self.path, k + '.pdb'))
                st = set()
                readPDB = open(os.path.join(self.path, k + '.pdb'), "rU")
                for record in SeqIO.parse(readPDB, "pdb-seqres"):
                    st.add(record.seq)
                seq.append(['>' + k + '\n' + ''.join([str(i) for i in st])])
            # add target seq
            t_seq_file = os.path.join(self.path, self.s.group(1) + '.pir')
            if os.path.exists(t_seq_file):
                t_seq = AlignIO.read(t_seq_file, 'pir')[0]
                seq.append(['>' + self.s.group(1) + '\n' + ''.join([i for i in t_seq.seq])])
            # clustal omega params
            params_omega = {}
            params_omega['title'] = ''
            # email
            params_omega['email'] = 'pymodel.v1@gmail.com'
            # dealign
            params_omega['dealign'] = 'False'
            # mbed clustering
            params_omega['mbediteration'] = 'True'
            # mbed guid tree
            params_omega['mbed'] = 'True'
            # output distance matrix
            params_omega['dismatout'] = 'True'
            # gtiterations
            params_omega['gtiterations'] = -1
            # hmmiterations
            params_omega['hmmiterations'] = -1
            # outfmt
            params_omega['outfmt'] = 'nexus'
            # seq type
            params_omega['stype'] = 'protein'
            # sequence
            params_omega['sequence'] = '\n'.join(i[0] for i in seq)
            # run clustalo
            run_clustalo.runClustalO.runclaustalo(params=params_omega)
            run_clustalo.runClustalO.getResult(self.s.group(1))
            self.msgLabel.setText('')
            self.compareBut.setEnabled(True)
            self.viewBut.setVisible(True)
        except Exception as er:
            self.msgLabel.setStyleSheet('color: red')
            self.msgLabel.setText('Error in comparing templates')
            self.compareBut.setEnabled(True)
            QApplication.processEvents()
            pass

        # MODELLER compare
        # file = open('info')
        # infoFile = json.load(file)
        # modeller_path = infoFile['MODELLER_PATH']
        # sys.path.insert(0, modeller_path)

        # for file in self.pdbL:
        #     if file[0] in self.pdb_index and os.path.exists(os.path.join(self.path, file[1] + '.pdb')):
        #         s = p.get_structure(os.path.join(self.path, file[1]), os.path.join(self.path, file[1] + '.pdb'))
        #         pdbCodeResolution.append([os.path.join(self.path, file[1] + '.pdb'), [chain.id for chain in s[0]]])

        # from modeller import log, environ, alignment, model, ModellerError
        #
        # #sys.stdout = open('log.txt', 'w')
        # log.verbose()
        # env = environ()
        # aln = alignment(env)
        # try:
        #     for i in pdbCodeResolution_unique:
        #         m = model(env, file=i[0], model_segment=('FIRST:' + i[1][0], 'LAST:' + i[1][len(i[1]) - 1]))
        #         aln.append_model(m, atom_files=i[0],
        #                          align_codes=os.path.basename(os.path.splitext(i[0])[0]) + ':' + ''.join(
        #                              [i for i in i[1]]))
        #     aln.malign()
        #     aln.malign3d()
        #     aln.compare_structures()
        #     aln.id_table(matrix_file='family.mat')
        #     env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)
        #     self.compare_window.show()
            #sys.stdout.close()
        # except ModellerError as er:
        #     print(er)
        #     pass

    def view_compare_result(self):
        # draw plot
        try:
            dnd_file = [i for i in os.listdir(self.path) if i.endswith('.dnd')]
            if len(dnd_file) > 0:
                for k, v in self.c_id.items():
                    with fileinput.FileInput(os.path.join(self.path, dnd_file[0]), inplace=True) as f:
                            for line in f:
                                print(line.replace(k, k + '_' + v))
                fig, ax = plt.subplots()
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.set_xlabel('Diversity')
                ax.set_title('Phylogenetic Tree')
                tree = Phylo.read(os.path.join(self.path, dnd_file[0]), 'newick')
                Phylo.draw(tree, axes=ax)
                plt.show()
        except Exception as er:
            self.viewBut.setText('View comparing result (Error)')
            print('Error:', er, 'Line {}.'.format(sys.exc_info()[-1].tb_lineno))
            pass


# def main():
# #     app = QApplication(sys.argv)
# #     form = SelectTemp()
# #     form.show()
# #     app.exec_()
# #
# #
# # if __name__ == '__main__':
# #     main()