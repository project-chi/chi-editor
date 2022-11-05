from PyQt6.QtGui import QAction, QIcon


class EraserButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Eraser")
        self.setIcon(QIcon("..\\resources\\eraser.png"))
