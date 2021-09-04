import os
import sys
import re
import json
import numpy as np
import matplotlib as plt
from matplotlib import cm
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from Bio import SeqIO, AlignIO
import tpl_view_align


class PlainTextEditor(QPlainTextEdit):
    def __init__(self, parent=None, **kwargs):
        super(PlainTextEditor, self).__init__(parent, **kwargs)
        font = QFont("Monospace")
        font.setPointSize(10)
        font.setStyleHint(QFont.TypeWriter)
        self.setFont(font)

    def zoom(self, delta):
        if delta < 0:
            self.zoomOut(1)
        elif delta > 0:
            self.zoomIn(1)

    def contextMenuEvent(self, event):

        menu = QMenu(self)
        copyAction = menu.addAction("Copy")
        clearAction = menu.addAction("Clear")
        zoominAction = menu.addAction("Zoom In")
        zoomoutAction = menu.addAction("Zoom Out")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == copyAction:
            self.copy()
        elif action == clearAction:
            self.clear()
        elif action == zoominAction:
            self.zoom(1)
        elif action == zoomoutAction:
            self.zoom(-1)


class AlignmentWidget(QWidget):
    def __init__(self, parent=None):
        super(AlignmentWidget, self).__init__(parent)
        l = QHBoxLayout(self)
        self.setLayout(l)
        self.m = QSplitter(self)
        l.addWidget(self.m)
        self.left = PlainTextEditor(self.m, readOnly=True)
        self.right = PlainTextEditor(self.m, readOnly=True)
        self.left.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.right.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.m.setSizes([200, 300])
        self.m.setStretchFactor(1, 2)


class SequencesViewer(QMainWindow, tpl_view_align.Ui_MainWindow):
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']
    file = ''
    hl = ''
    recs = None
    aln = None

    def __init__(self, parent=None):
        super(SequencesViewer, self).__init__(parent)
        self.setupUi(self)
        self.main = self.centralwidget
        self.setWindowTitle('Alignment Viewer')
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setGeometry(QtCore.QRect(200, 200, 960, 600))
        self.setMinimumHeight(150)
        self.zoomIn.clicked.connect(self.zoom_in)
        self.zoomOut.clicked.connect(self.zoom_out)
        self.alnFile.addItems(['Select'] + [i for i in os.listdir(self.path) if i.endswith('.aln')])
        self.alnFile.activated.connect(self.load_aln)
        self.browsAln.released.connect(self.brows_aln)
        self.start.textChanged.connect(self.show_alignment)
        self.end.textChanged.connect(self.show_alignment)
        self.search.textChanged.connect(self.search_seq)
        self.ed = PlainTextEditor(self, readOnly=True)
        self.ed.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.tabs.addTab(self.ed, 'Sequence')
        self.alnview = AlignmentWidget(self.main)
        self.tabs.addTab(self.alnview, 'Alignment')

    def scroll_top(self):
        vScrollBar = self.ed.verticalScrollBar()
        vScrollBar.triggerAction(QScrollBar.SliderToMinimum)

    def zoom_out(self):
        self.ed.zoom(-1)
        self.alnview.left.zoom(-1)
        self.alnview.right.zoom(-1)

    def zoom_in(self):
        self.ed.zoom(1)
        self.alnview.left.zoom(1)
        self.alnview.right.zoom(1)

    def brows_aln(self):

        file, _ = QFileDialog.getOpenFileName(self, 'Open File', './', filter="Alignment Files(*.aln);;All Files(*.*)")
        if not file:
            return
        try:
            self.recs = list(SeqIO.parse(file, 'pir'))
            self.show_seq()
            self.file = file
            self.scroll_top()
            self.alnview.left.clear()
            self.alnview.right.clear()
            self.start.setText('0')
            self.end.clear()
            self.end.setText(str(len(AlignIO.read(file, 'pir')[0])))
            self.show_alignment()
        except:
            pass

    def load_aln(self):
        if self.alnFile.currentText() != 'Select':
            try:
                alnFile = os.path.join(self.path, self.alnFile.currentText())
                self.recs = list(SeqIO.parse(alnFile, 'pir'))
                self.show_seq()
                self.file = alnFile
                self.scroll_top()
                self.alnview.left.clear()
                self.alnview.right.clear()
                self.start.setText('0')
                self.end.clear()
                self.end.setText(str(len(AlignIO.read(alnFile, 'pir')[0])))
                self.show_alignment()
            except:
                pass

    def search_seq(self):
        query = self.search.text().upper()
        aln = AlignIO.read(self.file, 'pir')
        if len(query.strip()) > 0:
            for i in aln:
                if query in i.seq:
                    self.hl = i.id
                    start = i.seq.index(query)
                    self.start.setText(str(start))
                    self.end.setText(str(start + len(query)))
                    break
                else:
                    self.hl = ''
                    self.start.setText('')
                    self.end.setText('')
                    continue
        else:
            self.hl = ''
            self.start.setText('0')
            self.end.setText(str(len(aln[0])))
        self.show_alignment()

    def show_seq(self):
        if self.recs is None:
            return
        self.ed.clear()
        for rec in self.recs:
            s = rec.format('pir')
            self.ed.insertPlainText(s)
        self.scroll_top()

    def show_alignment(self):
        try:
            self.alnview.left.clear()
            self.alnview.right.clear()
            if not re.search(r'\d*', self.start.text().strip()).group():
                start = 0
            else:
                start = int(self.start.text())
            if not re.search(r'\d*', self.end.text().strip()).group():
                end = 0
            else:
                end = int(self.end.text())
            colors = self.get_protein_colors()
            format = QtGui.QTextCharFormat()
            format.setBackground(QtGui.QBrush(QtGui.QColor('white')))
            cursorR = self.alnview.right.textCursor()
            cursorL = self.alnview.left.textCursor()
            aln = AlignIO.read(self.file, 'pir')
            lbls = np.arange(start, end + 1, 10)
            head = ''.join([('%-10s' % i) for i in lbls])
            cursorR.insertText(head)
            self.alnview.right.insertPlainText('\n')
            self.alnview.left.appendPlainText(' \n')
            for a in aln:
                seq = a.seq[start:end]
                if a.id == self.hl:
                    cursorL.insertHtml('<span style="background-color:#d6d30b;">%s</span>' % a.id)
                    self.alnview.left.insertPlainText('\n')
                else:
                    cursorL.insertHtml(a.id)
                    self.alnview.left.insertPlainText('\n')
                line = ''
                for aa in seq:
                   if aa in colors.keys():
                        c = colors[aa]
                        line += '<span style="background-color:%s;">%s</span>' % (c, aa)
                cursorR.insertHtml(line)
                self.alnview.right.insertPlainText('\n')
        except:
            pass

    @staticmethod
    def get_protein_colors(palette='tab20'):
        from Bio.PDB.Polypeptide import aa1
        aa1 = list(aa1)
        aa1.append('-')
        aa1.append('X')
        pal = cm.get_cmap(palette, 256)
        pal = [plt.colors.to_hex(i) for i in pal(np.linspace(0, 1, 20))]
        pal.append('white')
        pal.append('white')
        return {i: j for i, j in zip(aa1, pal)}


# def main():
#     app = QApplication(sys.argv)
#     sv = SequencesViewer()
#     sv.show()
#     app.exec_()
#
#
# if __name__ == '__main__':
#     main()
