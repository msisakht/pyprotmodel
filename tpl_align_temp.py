# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'tpl_align_temp.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(987, 577)
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(Form)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(10, 210, 971, 351))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox = QtWidgets.QGroupBox(self.verticalLayoutWidget_2)
        self.groupBox.setObjectName("groupBox")
        self.horizontalLayoutWidget_90 = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget_90.setGeometry(QtCore.QRect(0, 300, 961, 22))
        self.horizontalLayoutWidget_90.setObjectName("horizontalLayoutWidget_90")
        self.horizontalLayout_14 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_90)
        self.horizontalLayout_14.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_14.addItem(spacerItem)
        self.msgLabel2 = QtWidgets.QLabel(self.horizontalLayoutWidget_90)
        self.msgLabel2.setText("")
        self.msgLabel2.setObjectName("msgLabel2")
        self.horizontalLayout_14.addWidget(self.msgLabel2)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_14.addItem(spacerItem1)
        self.horizontalLayoutWidget_22 = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget_22.setGeometry(QtCore.QRect(10, 20, 171, 161))
        self.horizontalLayoutWidget_22.setObjectName("horizontalLayoutWidget_22")
        self.horizontalLayout_15 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_22)
        self.horizontalLayout_15.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_15.setObjectName("horizontalLayout_15")
        self.alignPDB = QtWidgets.QListWidget(self.horizontalLayoutWidget_22)
        self.alignPDB.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.alignPDB.setObjectName("alignPDB")
        self.horizontalLayout_15.addWidget(self.alignPDB)
        self.horizontalLayoutWidget_24 = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget_24.setGeometry(QtCore.QRect(310, 20, 191, 161))
        self.horizontalLayoutWidget_24.setObjectName("horizontalLayoutWidget_24")
        self.horizontalLayout_17 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_24)
        self.horizontalLayout_17.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_17.setObjectName("horizontalLayout_17")
        self.alignpdbChainAdded = QtWidgets.QListWidget(self.horizontalLayoutWidget_24)
        self.alignpdbChainAdded.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.alignpdbChainAdded.setObjectName("alignpdbChainAdded")
        self.horizontalLayout_17.addWidget(self.alignpdbChainAdded)
        self.horizontalLayoutWidget_23 = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget_23.setGeometry(QtCore.QRect(190, 20, 75, 161))
        self.horizontalLayoutWidget_23.setObjectName("horizontalLayoutWidget_23")
        self.horizontalLayout_16 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_23)
        self.horizontalLayout_16.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_16.setObjectName("horizontalLayout_16")
        self.alignChainPDB = QtWidgets.QListWidget(self.horizontalLayoutWidget_23)
        self.alignChainPDB.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.alignChainPDB.setObjectName("alignChainPDB")
        self.horizontalLayout_16.addWidget(self.alignChainPDB)
        self.horizontalLayoutWidget_89 = QtWidgets.QWidget(self.groupBox)
        self.horizontalLayoutWidget_89.setGeometry(QtCore.QRect(0, 320, 961, 31))
        self.horizontalLayoutWidget_89.setObjectName("horizontalLayoutWidget_89")
        self.horizontalLayout_114 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_89)
        self.horizontalLayout_114.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_114.setObjectName("horizontalLayout_114")
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_114.addItem(spacerItem2)
        self.alignBut = QtWidgets.QPushButton(self.horizontalLayoutWidget_89)
        self.alignBut.setObjectName("alignBut")
        self.horizontalLayout_114.addWidget(self.alignBut)
        self.viewAlign = QtWidgets.QPushButton(self.horizontalLayoutWidget_89)
        self.viewAlign.setEnabled(True)
        self.viewAlign.setObjectName("viewAlign")
        self.horizontalLayout_114.addWidget(self.viewAlign)
        self.ali_remove_file = QtWidgets.QPushButton(self.horizontalLayoutWidget_89)
        self.ali_remove_file.setObjectName("ali_remove_file")
        self.horizontalLayout_114.addWidget(self.ali_remove_file)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_114.addItem(spacerItem3)
        self.verticalLayoutWidget_3 = QtWidgets.QWidget(self.groupBox)
        self.verticalLayoutWidget_3.setGeometry(QtCore.QRect(266, 20, 41, 161))
        self.verticalLayoutWidget_3.setObjectName("verticalLayoutWidget_3")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_3)
        self.verticalLayout_5.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        spacerItem4 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_5.addItem(spacerItem4)
        self.alignAddBut = QtWidgets.QPushButton(self.verticalLayoutWidget_3)
        self.alignAddBut.setObjectName("alignAddBut")
        self.verticalLayout_5.addWidget(self.alignAddBut)
        self.addAll = QtWidgets.QPushButton(self.verticalLayoutWidget_3)
        self.addAll.setObjectName("addAll")
        self.verticalLayout_5.addWidget(self.addAll)
        self.alignMinBut = QtWidgets.QPushButton(self.verticalLayoutWidget_3)
        self.alignMinBut.setObjectName("alignMinBut")
        self.verticalLayout_5.addWidget(self.alignMinBut)
        spacerItem5 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_5.addItem(spacerItem5)
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(520, 20, 171, 151))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.localGlobal = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.localGlobal.setObjectName("localGlobal")
        self.gridLayout.addWidget(self.localGlobal, 1, 1, 1, 1)
        self.matrix = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.matrix.setObjectName("matrix")
        self.gridLayout.addWidget(self.matrix, 2, 1, 1, 1)
        self.alignType = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.alignType.setObjectName("alignType")
        self.gridLayout.addWidget(self.alignType, 0, 1, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_17.setObjectName("label_17")
        self.gridLayout.addWidget(self.label_17, 4, 0, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_15.setObjectName("label_15")
        self.gridLayout.addWidget(self.label_15, 2, 0, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_8.setObjectName("label_8")
        self.gridLayout.addWidget(self.label_8, 1, 0, 1, 1)
        self.gapOpen1d = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.gapOpen1d.setObjectName("gapOpen1d")
        self.gridLayout.addWidget(self.gapOpen1d, 4, 1, 1, 1)
        self.label_16 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_16.setObjectName("label_16")
        self.gridLayout.addWidget(self.label_16, 3, 0, 1, 1)
        self.gapFunction = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.gapFunction.setObjectName("gapFunction")
        self.gridLayout.addWidget(self.gapFunction, 3, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 0, 0, 1, 1)
        self.label_18 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_18.setObjectName("label_18")
        self.gridLayout.addWidget(self.label_18, 5, 0, 1, 1)
        self.gapExten1d = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.gapExten1d.setObjectName("gapExten1d")
        self.gridLayout.addWidget(self.gapExten1d, 5, 1, 1, 1)
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.groupBox)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(710, 20, 251, 151))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_20 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_20.setObjectName("label_20")
        self.gridLayout_2.addWidget(self.label_20, 0, 0, 1, 1)
        self.gapPenalt2dStright = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.gapPenalt2dStright.setObjectName("gapPenalt2dStright")
        self.gridLayout_2.addWidget(self.gapPenalt2dStright, 3, 1, 1, 1)
        self.gapPenalt2dCACA = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.gapPenalt2dCACA.setObjectName("gapPenalt2dCACA")
        self.gridLayout_2.addWidget(self.gapPenalt2dCACA, 4, 1, 1, 1)
        self.gapPenalt2dAccess = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.gapPenalt2dAccess.setObjectName("gapPenalt2dAccess")
        self.gridLayout_2.addWidget(self.gapPenalt2dAccess, 2, 1, 1, 1)
        self.label_25 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_25.setObjectName("label_25")
        self.gridLayout_2.addWidget(self.label_25, 4, 0, 1, 1)
        self.label_23 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_23.setObjectName("label_23")
        self.gridLayout_2.addWidget(self.label_23, 3, 0, 1, 1)
        self.gapPenalt2dHelix = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.gapPenalt2dHelix.setObjectName("gapPenalt2dHelix")
        self.gridLayout_2.addWidget(self.gapPenalt2dHelix, 0, 1, 1, 1)
        self.label_21 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_21.setObjectName("label_21")
        self.gridLayout_2.addWidget(self.label_21, 1, 0, 1, 1)
        self.gapPenalt2dBeta = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.gapPenalt2dBeta.setObjectName("gapPenalt2dBeta")
        self.gridLayout_2.addWidget(self.gapPenalt2dBeta, 1, 1, 1, 1)
        self.label_22 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_22.setObjectName("label_22")
        self.gridLayout_2.addWidget(self.label_22, 2, 0, 1, 1)
        self.label_26 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_26.setObjectName("label_26")
        self.gridLayout_2.addWidget(self.label_26, 5, 0, 1, 1)
        self.gapPenalt2dDstMin = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.gapPenalt2dDstMin.setObjectName("gapPenalt2dDstMin")
        self.gridLayout_2.addWidget(self.gapPenalt2dDstMin, 5, 1, 1, 1)
        self.gridLayoutWidget_3 = QtWidgets.QWidget(self.groupBox)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(10, 200, 951, 101))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.gridLayout_3.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_3.setObjectName("gridLayout_3")
        spacerItem6 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem6, 3, 2, 1, 1)
        spacerItem7 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem7, 3, 5, 1, 1)
        spacerItem8 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem8, 3, 8, 1, 1)
        self.similarityFlag = QtWidgets.QComboBox(self.gridLayoutWidget_3)
        self.similarityFlag.setObjectName("similarityFlag")
        self.gridLayout_3.addWidget(self.similarityFlag, 3, 7, 1, 1)
        self.matrixOffsetlabel = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.matrixOffsetlabel.setEnabled(False)
        self.matrixOffsetlabel.setObjectName("matrixOffsetlabel")
        self.gridLayout_3.addWidget(self.matrixOffsetlabel, 1, 9, 1, 1)
        self.rmsCutOff = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.rmsCutOff.setObjectName("rmsCutOff")
        self.gridLayout_3.addWidget(self.rmsCutOff, 0, 10, 1, 1)
        self.userMatrix = QtWidgets.QPushButton(self.gridLayoutWidget_3)
        self.userMatrix.setObjectName("userMatrix")
        self.gridLayout_3.addWidget(self.userMatrix, 2, 10, 1, 1)
        self.label_46 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_46.setObjectName("label_46")
        self.gridLayout_3.addWidget(self.label_46, 2, 9, 1, 1)
        self.matrixOffset = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.matrixOffset.setEnabled(False)
        self.matrixOffset.setObjectName("matrixOffset")
        self.gridLayout_3.addWidget(self.matrixOffset, 1, 10, 1, 1)
        self.overhangFactor = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.overhangFactor.setEnabled(False)
        self.overhangFactor.setObjectName("overhangFactor")
        self.gridLayout_3.addWidget(self.overhangFactor, 2, 7, 1, 1)
        self.label_37 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_37.setObjectName("label_37")
        self.gridLayout_3.addWidget(self.label_37, 3, 6, 1, 1)
        self.overHangFaclabel = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.overHangFaclabel.setEnabled(False)
        self.overHangFaclabel.setObjectName("overHangFaclabel")
        self.gridLayout_3.addWidget(self.overHangFaclabel, 2, 6, 1, 1)
        self.label_34 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_34.setObjectName("label_34")
        self.gridLayout_3.addWidget(self.label_34, 0, 3, 1, 1)
        self.gapResScore = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.gapResScore.setObjectName("gapResScore")
        self.gridLayout_3.addWidget(self.gapResScore, 1, 4, 1, 1)
        self.gapGapScore = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.gapGapScore.setObjectName("gapGapScore")
        self.gridLayout_3.addWidget(self.gapGapScore, 2, 4, 1, 1)
        self.autoOverhangLim = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.autoOverhangLim.setEnabled(False)
        self.autoOverhangLim.setObjectName("autoOverhangLim")
        self.gridLayout_3.addWidget(self.autoOverhangLim, 1, 7, 1, 1)
        self.overhang = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.overhang.setObjectName("overhang")
        self.gridLayout_3.addWidget(self.overhang, 3, 4, 1, 1)
        self.dendrogram = QtWidgets.QComboBox(self.gridLayoutWidget_3)
        self.dendrogram.setObjectName("dendrogram")
        self.gridLayout_3.addWidget(self.dendrogram, 2, 13, 1, 1)
        self.maxgapLen = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.maxgapLen.setObjectName("maxgapLen")
        self.gridLayout_3.addWidget(self.maxgapLen, 0, 4, 1, 1)
        self.useFit = QtWidgets.QComboBox(self.gridLayoutWidget_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.useFit.sizePolicy().hasHeightForWidth())
        self.useFit.setSizePolicy(sizePolicy)
        self.useFit.setObjectName("useFit")
        self.gridLayout_3.addWidget(self.useFit, 0, 13, 1, 1)
        self.label_33 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_33.setObjectName("label_33")
        self.gridLayout_3.addWidget(self.label_33, 3, 3, 1, 1)
        self.gapOpen3d = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gapOpen3d.sizePolicy().hasHeightForWidth())
        self.gapOpen3d.setSizePolicy(sizePolicy)
        self.gapOpen3d.setObjectName("gapOpen3d")
        self.gridLayout_3.addWidget(self.gapOpen3d, 2, 1, 1, 1)
        self.autoOverHangLimlabel = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.autoOverHangLimlabel.setEnabled(False)
        self.autoOverHangLimlabel.setObjectName("autoOverHangLimlabel")
        self.gridLayout_3.addWidget(self.autoOverHangLimlabel, 1, 6, 1, 1)
        self.label_48 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_48.setObjectName("label_48")
        self.gridLayout_3.addWidget(self.label_48, 2, 12, 1, 1)
        self.userMatrixType = QtWidgets.QComboBox(self.gridLayoutWidget_3)
        self.userMatrixType.setEnabled(False)
        self.userMatrixType.setObjectName("userMatrixType")
        self.gridLayout_3.addWidget(self.userMatrixType, 3, 10, 1, 1)
        self.label_51 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_51.setObjectName("label_51")
        self.gridLayout_3.addWidget(self.label_51, 0, 12, 1, 1)
        self.label_45 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_45.setObjectName("label_45")
        self.gridLayout_3.addWidget(self.label_45, 0, 9, 1, 1)
        self.label_29 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_29.setObjectName("label_29")
        self.gridLayout_3.addWidget(self.label_29, 1, 0, 1, 1)
        self.writeFit = QtWidgets.QComboBox(self.gridLayoutWidget_3)
        self.writeFit.setObjectName("writeFit")
        self.gridLayout_3.addWidget(self.writeFit, 1, 13, 1, 1)
        self.label_39 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_39.setObjectName("label_39")
        self.gridLayout_3.addWidget(self.label_39, 0, 6, 1, 1)
        self.label_28 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_28.setObjectName("label_28")
        self.gridLayout_3.addWidget(self.label_28, 0, 0, 1, 1)
        self.gapExten3d = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gapExten3d.sizePolicy().hasHeightForWidth())
        self.gapExten3d.setSizePolicy(sizePolicy)
        self.gapExten3d.setObjectName("gapExten3d")
        self.gridLayout_3.addWidget(self.gapExten3d, 3, 1, 1, 1)
        self.label_50 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_50.setObjectName("label_50")
        self.gridLayout_3.addWidget(self.label_50, 1, 12, 1, 1)
        self.label_35 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_35.setObjectName("label_35")
        self.gridLayout_3.addWidget(self.label_35, 1, 3, 1, 1)
        self.userMatrixTypelabel = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.userMatrixTypelabel.setEnabled(False)
        self.userMatrixTypelabel.setObjectName("userMatrixTypelabel")
        self.gridLayout_3.addWidget(self.userMatrixTypelabel, 3, 9, 1, 1)
        self.label_31 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_31.setObjectName("label_31")
        self.gridLayout_3.addWidget(self.label_31, 2, 0, 1, 1)
        self.autoOverhang = QtWidgets.QComboBox(self.gridLayoutWidget_3)
        self.autoOverhang.setObjectName("autoOverhang")
        self.gridLayout_3.addWidget(self.autoOverhang, 0, 7, 1, 1)
        self.label_30 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_30.setObjectName("label_30")
        self.gridLayout_3.addWidget(self.label_30, 3, 0, 1, 1)
        self.label_32 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_32.setObjectName("label_32")
        self.gridLayout_3.addWidget(self.label_32, 2, 3, 1, 1)
        self.gapPenalt2dDstPower = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gapPenalt2dDstPower.sizePolicy().hasHeightForWidth())
        self.gapPenalt2dDstPower.setSizePolicy(sizePolicy)
        self.gapPenalt2dDstPower.setObjectName("gapPenalt2dDstPower")
        self.gridLayout_3.addWidget(self.gapPenalt2dDstPower, 0, 1, 1, 1)
        self.gapPenalt2dStrucProfile = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.gapPenalt2dStrucProfile.sizePolicy().hasHeightForWidth())
        self.gapPenalt2dStrucProfile.setSizePolicy(sizePolicy)
        self.gapPenalt2dStrucProfile.setObjectName("gapPenalt2dStrucProfile")
        self.gridLayout_3.addWidget(self.gapPenalt2dStrucProfile, 1, 1, 1, 1)
        spacerItem9 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem9, 2, 11, 1, 1)
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(Form)
        self.groupBox_2.setGeometry(QtCore.QRect(10, 10, 971, 181))
        self.groupBox_2.setObjectName("groupBox_2")
        self.horizontalLayoutWidget_92 = QtWidgets.QWidget(self.groupBox_2)
        self.horizontalLayoutWidget_92.setGeometry(QtCore.QRect(0, 150, 971, 31))
        self.horizontalLayoutWidget_92.setObjectName("horizontalLayoutWidget_92")
        self.horizontalLayout_116 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_92)
        self.horizontalLayout_116.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_116.setObjectName("horizontalLayout_116")
        spacerItem10 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_116.addItem(spacerItem10)
        self.gen_ali = QtWidgets.QPushButton(self.horizontalLayoutWidget_92)
        self.gen_ali.setObjectName("gen_ali")
        self.horizontalLayout_116.addWidget(self.gen_ali)
        spacerItem11 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_116.addItem(spacerItem11)
        self.horizontalLayoutWidget_91 = QtWidgets.QWidget(self.groupBox_2)
        self.horizontalLayoutWidget_91.setGeometry(QtCore.QRect(0, 137, 971, 22))
        self.horizontalLayoutWidget_91.setObjectName("horizontalLayoutWidget_91")
        self.horizontalLayout_19 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_91)
        self.horizontalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_19.setObjectName("horizontalLayout_19")
        spacerItem12 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_19.addItem(spacerItem12)
        self.msg_gen = QtWidgets.QLabel(self.horizontalLayoutWidget_91)
        self.msg_gen.setText("")
        self.msg_gen.setObjectName("msg_gen")
        self.horizontalLayout_19.addWidget(self.msg_gen)
        spacerItem13 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_19.addItem(spacerItem13)
        self.gridLayoutWidget_5 = QtWidgets.QWidget(self.groupBox_2)
        self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 20, 661, 21))
        self.gridLayoutWidget_5.setObjectName("gridLayoutWidget_5")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.gridLayoutWidget_5)
        self.gridLayout_5.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.label = QtWidgets.QLabel(self.gridLayoutWidget_5)
        self.label.setObjectName("label")
        self.gridLayout_5.addWidget(self.label, 0, 0, 1, 1)
        self.gridLayoutWidget_6 = QtWidgets.QWidget(self.groupBox_2)
        self.gridLayoutWidget_6.setGeometry(QtCore.QRect(10, 40, 521, 98))
        self.gridLayoutWidget_6.setObjectName("gridLayoutWidget_6")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.gridLayoutWidget_6)
        self.gridLayout_6.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.seq_ali = QtWidgets.QTextEdit(self.gridLayoutWidget_6)
        self.seq_ali.setObjectName("seq_ali")
        self.gridLayout_6.addWidget(self.seq_ali, 0, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_2.setObjectName("label_2")
        self.gridLayout_6.addWidget(self.label_2, 0, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.gridLayoutWidget_6)
        self.label_3.setObjectName("label_3")
        self.gridLayout_6.addWidget(self.label_3, 1, 0, 1, 1)
        self.acc_ali = QtWidgets.QLineEdit(self.gridLayoutWidget_6)
        self.acc_ali.setObjectName("acc_ali")
        self.gridLayout_6.addWidget(self.acc_ali, 1, 1, 1, 1)
        self.gridLayoutWidget_7 = QtWidgets.QWidget(self.groupBox_2)
        self.gridLayoutWidget_7.setGeometry(QtCore.QRect(550, 40, 191, 51))
        self.gridLayoutWidget_7.setObjectName("gridLayoutWidget_7")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.gridLayoutWidget_7)
        self.gridLayout_7.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.label_4 = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.label_4.setObjectName("label_4")
        self.gridLayout_7.addWidget(self.label_4, 1, 0, 1, 1)
        self.label_49 = QtWidgets.QLabel(self.gridLayoutWidget_7)
        self.label_49.setObjectName("label_49")
        self.gridLayout_7.addWidget(self.label_49, 0, 0, 1, 1)
        self.seq = QtWidgets.QComboBox(self.gridLayoutWidget_7)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.seq.sizePolicy().hasHeightForWidth())
        self.seq.setSizePolicy(sizePolicy)
        self.seq.setObjectName("seq")
        self.gridLayout_7.addWidget(self.seq, 0, 1, 1, 1)
        self.browse_ali = QtWidgets.QPushButton(self.gridLayoutWidget_7)
        self.browse_ali.setObjectName("browse_ali")
        self.gridLayout_7.addWidget(self.browse_ali, 1, 1, 1, 1)
        self.gridLayoutWidget_4 = QtWidgets.QWidget(self.groupBox_2)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(750, 70, 211, 21))
        self.gridLayoutWidget_4.setObjectName("gridLayoutWidget_4")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.gridLayoutWidget_4)
        self.gridLayout_4.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.browse_label = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.browse_label.setText("")
        self.browse_label.setObjectName("browse_label")
        self.gridLayout_4.addWidget(self.browse_label, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.groupBox.setTitle(_translate("Form", "Align Templates"))
        self.alignBut.setText(_translate("Form", "Align"))
        self.viewAlign.setText(_translate("Form", "View Alignment"))
        self.ali_remove_file.setText(_translate("Form", "Remove"))
        self.alignAddBut.setText(_translate("Form", "Add"))
        self.addAll.setText(_translate("Form", "Add all"))
        self.alignMinBut.setText(_translate("Form", "Delete"))
        self.label_17.setToolTip(_translate("Form", "<html><head/><body><p>Gap creation penalty for sequence/sequence alignment.</p></body></html>"))
        self.label_17.setText(_translate("Form", "1D gap open:"))
        self.label_15.setText(_translate("Form", "Matrix:"))
        self.label_8.setToolTip(_translate("Form", "<html><head/><body><p>Whether to do local as opposed to global alignment.</p></body></html>"))
        self.label_8.setText(_translate("Form", "Local alignment:"))
        self.label_16.setToolTip(_translate("Form", "<html><head/><body><p>Whether or not to switch on functional gap penalty.</p></body></html>"))
        self.label_16.setText(_translate("Form", "Gap function:"))
        self.label_5.setToolTip(_translate("Form", "<html><head/><body><p>Pairwise, Tree and Progressive alignment.</p></body></html>"))
        self.label_5.setText(_translate("Form", "Alignment type:"))
        self.label_18.setToolTip(_translate("Form", "<html><head/><body><p>Gap extension penalty for sequence/sequence alignment.</p></body></html>"))
        self.label_18.setText(_translate("Form", "1D gap extension:"))
        self.label_20.setText(_translate("Form", "2D helix gap penalty:"))
        self.label_25.setText(_translate("Form", "2D CA-CA distance factor gap penalty:"))
        self.label_23.setText(_translate("Form", "2D Straightness gap penalty:"))
        self.label_21.setText(_translate("Form", "2D beta gap penalty:"))
        self.label_22.setText(_translate("Form", "2D Accessibility gap penalty:"))
        self.label_26.setText(_translate("Form", "2D DST min gap penalty:"))
        self.matrixOffsetlabel.setToolTip(_translate("Form", "<html><head/><body><p>Substitution matrix offset for local alignment.</p></body></html>"))
        self.matrixOffsetlabel.setText(_translate("Form", "Matrix offset:"))
        self.userMatrix.setText(_translate("Form", "Browse"))
        self.label_46.setText(_translate("Form", "User matrix:"))
        self.label_37.setToolTip(_translate("Form", "<html><head/><body><p>If yes, the SALIGN command does not convert numbers into a distance sense.</p></body></html>"))
        self.label_37.setText(_translate("Form", "Similarity flag:"))
        self.overHangFaclabel.setToolTip(_translate("Form", "<html><head/><body><p>Factor to multiply sequence length difference with when auto overhang is yes.</p></body></html>"))
        self.overHangFaclabel.setText(_translate("Form", "Overhang factor:"))
        self.label_34.setText(_translate("Form", "Max gap length:"))
        self.label_33.setToolTip(_translate("Form", "<html><head/><body><p>Un-penalized overhangs in protein comparisons.</p></body></html>"))
        self.label_33.setText(_translate("Form", "Overhang:"))
        self.autoOverHangLimlabel.setToolTip(_translate("Form", "<html><head/><body><p>Auto overhang effective if sequence length difference greater than this parameter.</p></body></html>"))
        self.autoOverHangLimlabel.setText(_translate("Form", "Auto overhang limit:"))
        self.label_48.setText(_translate("Form", "Dendrogram:"))
        self.label_51.setToolTip(_translate("Form", "<html><head/><body><p>Whether to do pairwise least-squares fitting or ALIGN2D alignment.</p></body></html>"))
        self.label_51.setText(_translate("Form", "Use fit:"))
        self.label_45.setToolTip(_translate("Form", "<html><head/><body><p>Cutoffs for RMS, DRMS, Alpha, Phi, Psi, Omega, chi1, chi2, chi3, chi4 and chi5.</p></body></html>"))
        self.label_45.setText(_translate("Form", "RMS cut off:"))
        self.label_29.setText(_translate("Form", "2D structure profile gap penalty:"))
        self.label_39.setToolTip(_translate("Form", "<html><head/><body><p>Overhang values made dependent on sequence length difference.</p></body></html>"))
        self.label_39.setText(_translate("Form", "Auto overhang:"))
        self.label_28.setText(_translate("Form", "2D DST power gap penalty:"))
        self.label_50.setText(_translate("Form", "Write fit file:"))
        self.label_35.setText(_translate("Form", "Gap-residue score:"))
        self.userMatrixTypelabel.setText(_translate("Form", "Matrix type:"))
        self.label_31.setText(_translate("Form", "3D gap open:"))
        self.label_30.setText(_translate("Form", "3D gap extension:"))
        self.label_32.setText(_translate("Form", "Gap-gap score:"))
        self.groupBox_2.setTitle(_translate("Form", "Sequence Preparation"))
        self.gen_ali.setText(_translate("Form", "Generate"))
        self.label.setText(_translate("Form", "If you didn\'t use the PyModel to search for templates, please enter your sequences or accession numbers to generate the PIR files."))
        self.label_2.setText(_translate("Form", "Sequence:"))
        self.label_3.setToolTip(_translate("Form", "<html><head/><body><p>UniProt, NCBI or PDB accession numbers.</p></body></html>"))
        self.label_3.setText(_translate("Form", "Accession number(s):"))
        self.label_4.setText(_translate("Form", "Or"))
        self.label_49.setText(_translate("Form", "Select sequence file:"))
        self.browse_ali.setText(_translate("Form", "Browse file"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
