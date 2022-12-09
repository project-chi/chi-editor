from math import fabs, atan2, sin
from pathlib import Path

from PyQt6.QtWidgets import QGraphicsItem, QGraphicsPixmapItem
from PyQt6.QtGui import QPixmap, QPainter
from PyQt6.QtCore import QPointF


class Line(QGraphicsPixmapItem):
    vertex1: QGraphicsItem
    vertex2: QGraphicsItem
    end1: QPointF
    end2: QPointF
    width: float
    height: float
    MAX_WIDTH = 30

    def __init__(self, *args, vertex1: QGraphicsItem, vertex2: QGraphicsItem, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.vertex1 = vertex1
        self.vertex2 = vertex2

        # centers of vertices
        self.end1 = QPointF(vertex1.x() + vertex1.boundingRect().width() / 2,
                            vertex1.y() + vertex1.boundingRect().height() / 2)
        self.end1 = QPointF(vertex2.x() + vertex2.boundingRect().width() / 2,
                            vertex2.y() + vertex2.boundingRect().height() / 2)

        # line will be rotated around its vertices
        self.setTransformOriginPoint(self.end1)

        self.setRotation(atan2(self.end1.y() - self.end2.y(), self.end1.x() - self.end2.x()))

        # min of MAX_WIDTH and boarders of vertex items
        self.width = min([self.MAX_WIDTH, vertex1.boundingRect().width(), vertex1.boundingRect().height(),
                            vertex2.boundingRect().width(), vertex2.boundingRect().height()])

        # hypotenuse of right triangle of vertices, rotation is an angle between them
        self.height = fabs(self.end1.y() - self.end1.y()) / sin(self.rotation())

        # height is just distance between y's
        self.setPixmap(QPixmap(self.width, self.height))

    def paint(self):
        painter = QPainter(self.pixmap())
        rect = self.texture.rect()
        painter.drawLine(rect.x() + rect.width() / 2, rect.y(),
                         rect.x() + rect.width() / 2, rect.y() + rect.height())
