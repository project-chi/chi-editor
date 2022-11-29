from PyQt6 import QtGui
from PyQt6.QtCore import QPoint
from PyQt6.QtWidgets import QLabel


class Canvas(QLabel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.canvas = QtGui.QPixmap(400, 300)
        self.canvas.fill(QtGui.QColor("white"))
        self.setPixmap(self.canvas)

    def mousePressEvent(self, ev: QtGui.QMouseEvent) -> None:
        text_widget = DragLabel(self)
        text_widget.setText("Some text")
        text_widget.move(ev.pos())
        text_widget.show()


class DragLabel(QLabel):
    mouse_move_start = False
    delta = QPoint()

    def __init__(self, parent):
        super().__init__(parent)

    def mousePressEvent(self, ev: QtGui.QMouseEvent) -> None:
        self.mouse_move_start = True
        self.delta = self.pos() - ev.scenePosition().toPoint()

    def mouseMoveEvent(self, ev: QtGui.QMouseEvent) -> None:
        if self.mouse_move_start:
            self.move(ev.scenePosition().toPoint() + self.delta)

    def mouseReleaseEvent(self, ev: QtGui.QMouseEvent) -> None:
        self.mouse_move_start = False
