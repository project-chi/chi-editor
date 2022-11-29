from PyQt6 import QtGui
from PyQt6.QtCore import QPoint
from PyQt6.QtWidgets import QLabel


class Canvas(QLabel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.canvas = QtGui.QPixmap(400, 300)
        self.canvas.fill(QtGui.QColor("white"))
        self.setPixmap(self.canvas)

    def mousePressEvent(self, ev: QtGui.QMouseEvent) -> None:
        text_widget = DragLabel(self, text='Some text')
        text_widget.move(ev.pos())
        text_widget.show()


class DragLabel(QLabel):
    is_moving: bool
    delta: QPoint
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_moving = False
        self.delta = QPoint()

    def mousePressEvent(self, ev: QtGui.QMouseEvent) -> None:
        self.is_moving = True
        self.delta = self.pos() - ev.scenePosition().toPoint()

    def mouseMoveEvent(self, ev: QtGui.QMouseEvent) -> None:
        if self.is_moving:
            self.move(ev.scenePosition().toPoint() + self.delta)

    def mouseReleaseEvent(self, ev: QtGui.QMouseEvent) -> None:
        self.is_moving = False
