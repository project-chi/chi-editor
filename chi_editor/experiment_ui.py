from PyQt6.QtCore import QRect, QCoreApplication, QMetaObject
from PyQt6.QtWidgets import QGraphicsView, QRadioButton, QButtonGroup


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(538, 269)
        self.graphicsView = QGraphicsView(Dialog)
        self.graphicsView.setGeometry(QRect(130, 10, 371, 221))
        self.graphicsView.setObjectName("graphicsView")
        self.retranslateUi(Dialog)
        QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
