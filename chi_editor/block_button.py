from PyQt6.QtGui import QAction, QIcon


class BlockButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Block")
        self.setIcon(QIcon("..\\resources\\block.png"))
