from PyQt6.QtGui import QAction, QIcon


class AtomButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Atom")
        self.setIcon(QIcon("..\\resources\\atom.png"))
