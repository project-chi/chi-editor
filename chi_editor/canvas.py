from PyQt6.QtGui import QColor, QPixmap, QMouseEvent, QPalette
from PyQt6.QtWidgets import QLabel, QLineEdit


class Canvas(QLabel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.canvas = QPixmap(400, 300)
        self.canvas.fill(QColor("white"))
        self.setPixmap(self.canvas)

    def mousePressEvent(self, ev: QMouseEvent) -> None:
        text_widget = QLabel("Some text", self)
        text_widget.move(ev.pos())
        text_widget.show()
