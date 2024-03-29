from math import atan2, cos, degrees, fabs, radians
from typing import TYPE_CHECKING

from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QColor, QPen, QPixmap
from PyQt6.QtWidgets import QGraphicsEllipseItem, QGraphicsItem, QGraphicsPixmapItem

if TYPE_CHECKING:
    from typing import ClassVar

    from PyQt6.QtCore import QRectF
    from PyQt6.QtGui import QPainter

    from chi_editor.bases.alpha_atom import AlphaAtom


class Line(QGraphicsPixmapItem):
    """Connects two graphics items with a line."""
    width: "ClassVar[float]" = 30.0

    vertex1: "AlphaAtom"
    vertex2: "QGraphicsItem"
    height: "float"
    multiplicity: "int"

    def __init__(
        self,
        start: "AlphaAtom",
        end: "QGraphicsItem | QPointF",
        *args, **kwargs,
    ) -> "None":
        super().__init__(*args, **kwargs)
        self.vertex1 = start

        if isinstance(end, QGraphicsItem):
            self.vertex2 = end
            self.vertex1.molecule.update_atoms()
        else:
            self.vertex2 = QGraphicsEllipseItem(0, 0, 0, 0)
            self.vertex2.setPos(end)

        self.setShapeMode(QGraphicsPixmapItem.ShapeMode.BoundingRectShape)

        self.setPos(
            self.vertex1.sceneBoundingRect().center() - QPointF(self.width / 2, 0)
        )

        self.update_pixmap(self.vertex2)

    def set_v2(self, end: "QGraphicsItem") -> "None":
        self.vertex2 = end

    def remove(self) -> "None":
        self.vertex1.lines.remove(self)
        self.vertex2.lines.remove(self)
        if self.scene():
            self.scene().removeItem(self)

    def update_pixmap(
        self,
        moved_vertex: "QGraphicsItem",
        following_mouse: "bool" = False,
    ) -> "None":
        moved_point = moved_vertex.sceneBoundingRect().center()
        self.setRotation(0)
        self.setScale(1)

        if moved_vertex is self.vertex2 or following_mouse:
            static_point = self.vertex1.sceneBoundingRect().center()
        else:
            static_point = self.vertex2.sceneBoundingRect().center()
        self.setTransformOriginPoint(self.mapFromScene(static_point))

        self.setPos(static_point - QPointF(self.width / 2, 0))

        self.setRotation(
            -1
            * degrees(
                atan2(
                    moved_point.x() - static_point.x(),
                    moved_point.y() - static_point.y(),
                )
            )
        )

        if self.rotation() % 360 == 90 or self.rotation() % 360 == 270:
            self.height = fabs(static_point.x() - moved_point.x())
        else:
            self.height = fabs(
                (static_point.y() - moved_point.y()) / cos(radians(self.rotation()))
            )

        self.setPixmap(QPixmap(int(self.width), int(self.height)))

        self.update()

    def paint(self, painter: "QPainter", *args) -> "None":
        painter.save()
        self.paint_line(painter)
        painter.restore()

    def paint_line(self, painter: "QPainter") -> "None":
        """This method determines how line is drawn.

        It should be overriden by children classes.
        """
        pen = QPen(QColor("black"), 3)
        painter.setPen(pen)

        painter.drawLine(
            QPointF(self.width / 2, 0), QPointF(self.width / 2, self.height)
        )

    def boundingRect(self) -> "QRectF":
        return super().boundingRect()
