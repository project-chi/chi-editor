from PyQt6.QtGui import QAction, QIcon


class StructureButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Structure")
        self.setIcon(QIcon("..\\resources\\structure.png"))
