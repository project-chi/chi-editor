from PyQt6.QtGui import QColor, QPixmap
from PyQt6.QtWidgets import QLabel


class Canvas(QLabel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.canvas = QPixmap(400, 300)
        self.canvas.fill(QColor("white"))
        self.setPixmap(self.canvas)
