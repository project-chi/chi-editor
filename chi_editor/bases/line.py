from math import fabs, atan2, cos, degrees, radians

from PyQt6.QtWidgets import QGraphicsItem, QGraphicsPixmapItem, QWidget, QStyleOptionGraphicsItem
from PyQt6.QtGui import QPixmap, QPainter
from PyQt6.QtCore import QPointF, QRectF


class Line(QGraphicsPixmapItem):
    vertex1: QGraphicsItem
    vertex2: QGraphicsItem
    width: float
    height: float
    MAX_WIDTH = 30

    def __init__(self, vertex1: QGraphicsItem, vertex2: QGraphicsItem, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.vertex1 = vertex1
        self.vertex2 = vertex2

        self.setShapeMode(QGraphicsPixmapItem.ShapeMode.BoundingRectShape)

        self.setPos(vertex1.sceneBoundingRect().center() - QPointF(vertex1.sceneBoundingRect().width() / 2, 0))

        # min of MAX_WIDTH and boarders of vertex items
        self.width = min([self.MAX_WIDTH, vertex1.boundingRect().width(), vertex1.boundingRect().height(),
                          vertex2.boundingRect().width(), vertex2.boundingRect().height()])

        self.__update_pixmap(self.vertex2)

    # recalculate height and rotation of pixmap
    def __update_pixmap(self, moved_vertex: QGraphicsItem) -> None:
        # simple deduction of static vertex
        moved_point = moved_vertex.sceneBoundingRect().center()

        if moved_vertex is self.vertex1:
            static_point = self.vertex2.sceneBoundingRect().center()
        else:
            static_point = self.vertex1.sceneBoundingRect().center()

        # line will be rotated around its vertices
        self.setTransformOriginPoint(self.mapFromScene(static_point))

        self.setRotation(-1 * degrees(atan2(moved_point.x() - static_point.x(), moved_point.y() - static_point.y())))

        # hypotenuse of right triangle of vertices, rotation is an angle between them
        if self.rotation() == 0:
            self.height = fabs(static_point.y() - moved_point.y())
        else:
            self.height = fabs((static_point.y() - moved_point.y()) / cos(radians(self.rotation())))

        # height is just distance between y's
        self.setPixmap(QPixmap(int(self.width), int(self.height)))

        self.update()

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem = None, widget: QWidget = None) -> None:
        pen = painter.pen()
        pen.setWidth(3)
        painter.setPen(pen)

        # draw straight line in the parallel to y-axis
        painter.drawLine(QPointF(self.width / 2, 0),
                         QPointF(self.width / 2, self.height))

    def boundingRect(self) -> QRectF:
        return super(Line, self).boundingRect()
