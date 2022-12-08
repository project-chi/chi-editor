from PyQt6.QtCore import QRect, QCoreApplication, QMetaObject
from PyQt6.QtWidgets import QGraphicsView, QRadioButton, QButtonGroup


class UIDialog:
    def setup_ui(self, dialog):
        dialog.setObjectName("Dialog")
        dialog.resize(538, 269)
        self.graphicsView = QGraphicsView(dialog)
        self.graphicsView.setGeometry(QRect(130, 10, 371, 221))
        self.graphicsView.setObjectName("graphicsView")
        self.retranslate_ui(dialog)
        QMetaObject.connectSlotsByName(dialog)

    def retranslate_ui(self, dialog):
        _translate = QCoreApplication.translate
        dialog.setWindowTitle(_translate("Dialog", "Dialog"))
