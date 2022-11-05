from PyQt6.QtGui import QAction, QIcon


class BondButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Bond")
        self.setIcon(QIcon("..\\resources\\bond.png"))
