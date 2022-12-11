from PyQt6.QtGui import QPen, QBrush, QColor
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem
from PyQt6.QtCore import QRectF

from ..canvas import Canvas
from ..constants import RESOURCES
from .line import Line


class AlphaAtom(QGraphicsItem):
    pen: QPen = QPen(QColor("black"), 0)
    brush: QBrush = QBrush(QColor("black"))
    rect: QRectF = QRectF(0, 0, 100, 100)
    lines: list[Line] = []

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setFlags(self.GraphicsItemFlag.ItemSendsScenePositionChanges)

    def addLine(self, newLine: Line):
        for existing in self.lines:
            if (existing.vertex1, existing.vertex2) == (newLine.vertex1, newLine.vertex2):
                # another line with the same control points already exists
                return False
        self.lines.append(newLine)
        return True

    def removeLine(self, line: Line):
        for existing in self.lines:
            if (existing.vertex1, existing.vertex2) == (line.vertex1, line.vertex2):
                self.scene().removeItem(existing)
                self.lines.remove(existing)
                return True
        return False

    def itemChange(self, change, value):
        for line in self.lines:
            line.update_pixmap(self)
        return super().itemChange(change, value)