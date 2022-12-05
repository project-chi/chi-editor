from PyQt6 import QtGui
from PyQt6.QtCore import QPoint
from PyQt6.QtGui import QBrush, QColor, QPen, QFont, QTextItem
from PyQt6.QtWidgets import QLabel, QLineEdit, QGraphicsView, QGraphicsScene, QGraphicsItem, QGraphicsSceneMouseEvent, \
    QGraphicsTextItem
from PyQt6.QtGui import QPen, QBrush, QColor, QFont


class Canvas(QGraphicsScene):
    shapes: list[QGraphicsItem]
    opt: str

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setSceneRect(-100, -100, 200, 200)
        self.shapes = []
        self.opt = ""

    def setOption(self, opt):
        self.opt = opt

    def toggleDragMode(self):
        if self.dragMode() == QtGui.QGraphicsView.ScrollHandDrag:
            self.setDragMode(QtGui.QGraphicsView.NoDrag)
        elif not self._photo.pixmap().isNull():
            self.setDragMode(QtGui.QGraphicsView.ScrollHandDrag)

    # def generateMode(self):
    #     if self.

    def mousePressEvent(self, event):
        if self.opt == "Generate":
            text = QGraphicsTextItem()
            text.setPos(event.scenePos())
            text.setPlainText("sadasdasd")
            text.setFont(QFont("Impact"))
            text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
            self.shapes.append(text)
            self.addItem(text)
        else:
            super(Canvas, self).mousePressEvent(event)

