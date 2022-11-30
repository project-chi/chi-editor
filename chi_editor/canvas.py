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


class Canvas(QLabel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.canvas = QtGui.QPixmap(400, 300)
        self.canvas.fill(QtGui.QColor("white"))
        self.setPixmap(self.canvas)

    def mousePressEvent(self, ev: QtGui.QMouseEvent) -> None:
        text_widget = DragLabel(self, text="Default text")
        text_widget.move(ev.pos())
        text_widget.show()

        text_line = QLineEdit(self)
        text_line.show()
        text_line.setFocus()
        text_line.move(ev.pos())
        text_line.returnPressed.connect(lambda: self.read_line(text_line, text_widget))

    def read_line(self, text_line: QLineEdit, text_widget: DragLabel) -> None:
        text_widget.setText(text_line.text())
        text_widget.adjustSize()

        text_line.close()
