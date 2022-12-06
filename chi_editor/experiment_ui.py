from PyQt6.QtCore import QRect, QCoreApplication, QMetaObject
from PyQt6.QtWidgets import QGraphicsView, QRadioButton, QButtonGroup


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(538, 269)
        self.graphicsView = QGraphicsView(Dialog)
        self.graphicsView.setGeometry(QRect(130, 10, 371, 221))
        self.graphicsView.setObjectName("graphicsView")
        self.buttonGroup = QButtonGroup
        self.radioButton = QRadioButton(Dialog)
        self.radioButton.setGeometry(QRect(20, 30, 82, 31))
        self.radioButton.setObjectName("radioButton")
        self.radioButton_2 = QRadioButton(Dialog)
        self.radioButton_2.setGeometry(QRect(20, 80, 82, 17))
        self.radioButton_2.setObjectName("radioButton_2")
        self.retranslateUi(Dialog)
        QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.radioButton.setText(_translate("Dialog", "Generate"))
        self.radioButton_2.setText(_translate("Dialog", "Select"))