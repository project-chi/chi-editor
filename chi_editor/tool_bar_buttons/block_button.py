from PyQt6.QtGui import QAction


class BlockButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Block")
