# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'tpl_view_align.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.tabs = QtWidgets.QTabWidget(self.centralwidget)
        self.tabs.setEnabled(True)
        self.tabs.setObjectName("tabs")
        self.horizontalLayout.addWidget(self.tabs)
        self.sidebarLayout = QtWidgets.QVBoxLayout()
        self.sidebarLayout.setObjectName("sidebarLayout")
        self.sidebar = QtWidgets.QWidget(self.centralwidget)
        self.sidebar.setMinimumSize(QtCore.QSize(180, 0))
        self.sidebar.setObjectName("sidebar")
        self.formLayoutWidget = QtWidgets.QWidget(self.sidebar)
        self.formLayoutWidget.setGeometry(QtCore.QRect(10, 120, 161, 25))
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.formLayout = QtWidgets.QFormLayout(self.formLayoutWidget)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.formLayout.setObjectName("formLayout")
        self.zoomIn = QtWidgets.QPushButton(self.formLayoutWidget)
        self.zoomIn.setObjectName("zoomIn")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.zoomIn)
        self.zoomOut = QtWidgets.QPushButton(self.formLayoutWidget)
        self.zoomOut.setObjectName("zoomOut")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.zoomOut)
        self.label = QtWidgets.QLabel(self.sidebar)
        self.label.setGeometry(QtCore.QRect(10, 0, 161, 16))
        self.label.setObjectName("label")
        self.alnFile = QtWidgets.QComboBox(self.sidebar)
        self.alnFile.setGeometry(QtCore.QRect(10, 20, 161, 22))
        self.alnFile.setObjectName("alnFile")
        self.browsAln = QtWidgets.QPushButton(self.sidebar)
        self.browsAln.setGeometry(QtCore.QRect(10, 50, 161, 23))
        self.browsAln.setObjectName("browsAln")
        self.horizontalLayoutWidget = QtWidgets.QWidget(self.sidebar)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 80, 160, 31))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.start = QtWidgets.QLineEdit(self.horizontalLayoutWidget)
        self.start.setText("")
        self.start.setObjectName("start")
        self.horizontalLayout_3.addWidget(self.start)
        self.end = QtWidgets.QLineEdit(self.horizontalLayoutWidget)
        self.end.setObjectName("end")
        self.horizontalLayout_3.addWidget(self.end)
        self.sidebarLayout.addWidget(self.sidebar)
        self.horizontalLayout.addLayout(self.sidebarLayout)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabs.setCurrentIndex(-1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.zoomIn.setText(_translate("MainWindow", "Zoom In"))
        self.zoomOut.setText(_translate("MainWindow", "Zoom Out"))
        self.label.setText(_translate("MainWindow", "Alignment file"))
        self.browsAln.setText(_translate("MainWindow", "Browse file"))
        self.start.setPlaceholderText(_translate("MainWindow", "Start"))
        self.end.setPlaceholderText(_translate("MainWindow", "End"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
