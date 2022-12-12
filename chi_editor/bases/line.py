from math import fabs, atan2, cos, degrees, radians

from PyQt6.QtWidgets import QGraphicsItem, QGraphicsPixmapItem, QWidget, QStyleOptionGraphicsItem, QGraphicsEllipseItem
from PyQt6.QtGui import QPixmap, QPainter, QPen, QColor
from PyQt6.QtCore import QPointF, QRectF


class Line(QGraphicsPixmapItem):
    """
    Class that connects two graphics items with a pixmap (presumably line).
    Parent class for all chemical bonds.

    Pixmap item is chosen instead of line item because chemical bond can be different (double, dotted, bold-wedged)
    """
    vertex1: QGraphicsItem
    vertex2: QGraphicsItem
    width: float
    height: float
    MAX_WIDTH = 30

    def __init__(self, start: QGraphicsItem, end: QGraphicsItem | QPointF, *args, **kwargs) \
            -> None:
        super().__init__(*args, **kwargs)
        self.vertex1 = start

        if isinstance(end, QGraphicsItem):
            self.vertex2 = end
            self.width = min([self.MAX_WIDTH, start.boundingRect().width(), start.boundingRect().height(),
                              end.boundingRect().width(), end.boundingRect().height()])
        else:
            self.vertex2 = QGraphicsEllipseItem(0, 0, 0, 0)
            self.vertex2.setPos(end)
            self.width = min([self.MAX_WIDTH, start.boundingRect().width(), start.boundingRect().height()])

        self.setShapeMode(QGraphicsPixmapItem.ShapeMode.BoundingRectShape)

        self.setPos(self.vertex1.sceneBoundingRect().center() - QPointF(self.width / 2, 0))

        self.update_pixmap(self.vertex2)

    def setV2(self, end: QGraphicsItem) -> None:
        self.vertex2 = end

    # recalculate height and rotation of pixmap
    def update_pixmap(self, moved_vertex: QGraphicsItem, following_mouse: bool = False) -> None:
        # simple deduction of static vertex
        moved_point = moved_vertex.sceneBoundingRect().center()
        self.setRotation(0)
        self.setScale(1)

        # moved_vertex is second vertex or bond follows the mouse
        if moved_vertex is self.vertex2 or following_mouse:
            static_point = self.vertex1.sceneBoundingRect().center()
        else:
            static_point = self.vertex2.sceneBoundingRect().center()
        self.setPos(static_point - QPointF(self.width / 2, 0))

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
        painter.save()
        self.paintLine(painter)
        painter.restore()

    def paintLine(self, painter: QPainter) -> None:
        """
        This method determines how line is drawn.
        It should be overrided by children classes
        """
        pen = QPen(QColor("black"), 3)
        painter.setPen(pen)

        # draw straight line
        painter.drawLine(QPointF(self.width / 2, 0),
                         QPointF(self.width / 2, self.height))

    def boundingRect(self) -> QRectF:
        return super(Line, self).boundingRect()
