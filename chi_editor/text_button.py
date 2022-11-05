from PyQt6.QtGui import QAction, QIcon


class TextButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Text")
        self.setIcon(QIcon("..\\resources\\text.png"))
