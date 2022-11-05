from PyQt6.QtGui import QAction


class AtomButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Atom")
