from PyQt6.QtGui import QAction, QIcon


class HandButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Hand")
        self.setIcon(QIcon("..\\resources\\hand.png"))
