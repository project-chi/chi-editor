from math import fabs, atan2, sin, degrees, radians
from pathlib import Path

from PyQt6.QtWidgets import QGraphicsItem, QGraphicsPixmapItem, QWidget, QStyleOptionGraphicsItem
from PyQt6.QtGui import QPixmap, QPainter
from PyQt6.QtCore import QPointF, QRectF


class Line(QGraphicsPixmapItem):
    vertex1: QGraphicsItem
    vertex2: QGraphicsItem
    width: float
    height: float
    MAX_WIDTH = 30

    def __init__(self, *args, vertex1: QGraphicsItem, vertex2: QGraphicsItem, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.vertex1 = vertex1
        self.vertex2 = vertex2

        self.setPos(vertex1.boundingRect().center())

        # min of MAX_WIDTH and boarders of vertex items
        self.width = min([self.MAX_WIDTH, vertex1.boundingRect().width(), vertex1.boundingRect().height(),
                          vertex2.boundingRect().width(), vertex2.boundingRect().height()])

        self.__update_pixmap(self.vertex2)

    # recalculate height and rotation of pixmap
    def __update_pixmap(self, moved_vertex: QGraphicsItem) -> None:
        # if (moved_vertex is not self.vertex1) or (moved_vertex is not self.vertex2):
        #     raise ValueError(f'Can\'t call __update_pixmap with {moved_vertex} argument')

        # simple deduction of static vertex
        moved_point = moved_vertex.boundingRect().center()
        if moved_vertex is self.vertex1:
            static_point = self.vertex2.boundingRect().center()
        else:
            static_point = self.vertex1.boundingRect().center()

        # line will be rotated around its vertices
        self.setTransformOriginPoint(static_point)

        self.setRotation(degrees(atan2(static_point.y() - moved_point.y(), static_point.x() - moved_point.x())))

        # hypotenuse of right triangle of vertices, rotation is an angle between them
        self.height = fabs((static_point.y() - moved_point.y()) / sin(radians(self.rotation())))

        # height is just distance between y's
        self.setPixmap(QPixmap(int(self.width), int(self.height)))

        self.update()

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem = None, widget: QWidget = None) -> None:
        rect = self.pixmap().rect()
        pen = painter.pen()
        pen.setWidth(5)
        painter.setPen(pen)

        # draw straight line in the parallel to y-axis
        painter.drawLine(QPointF(rect.x() + rect.width() / 2, rect.y()),
                         QPointF(rect.x() + rect.width() / 2, rect.y() + rect.height()))

    def boundingRect(self) -> QRectF:
        return QRectF(self.pixmap().rect())
