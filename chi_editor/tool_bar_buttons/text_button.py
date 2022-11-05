from PyQt6.QtGui import QAction


class TextButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Text")
