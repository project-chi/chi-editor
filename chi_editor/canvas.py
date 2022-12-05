from PyQt6 import QtGui
from PyQt6.QtCore import QPoint
from PyQt6.QtWidgets import QLabel, QLineEdit


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


class MyLineEdit(QLineEdit):
    def focusOutEvent(self, a0: QtGui.QFocusEvent) -> None:
        self.close()


class Canvas(QLabel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.canvas = QtGui.QPixmap(400, 300)
        self.canvas.fill(QtGui.QColor("white"))
        self.setPixmap(self.canvas)

    def mousePressEvent(self, ev: QtGui.QMouseEvent) -> None:
        text_line = MyLineEdit(self)
        text_line.show()
        text_line.setFocus()
        text_line.move(ev.pos())
        text_line.returnPressed.connect(lambda: self.read_line(text_line))

    def sth(self, event):
        print("lol")

    def read_line(self, text_line: QLineEdit) -> None:
        text_widget = DragLabel(self, text=text_line.text())
        text_widget.move(text_line.pos())
        text_widget.show()

        text_line.close()
