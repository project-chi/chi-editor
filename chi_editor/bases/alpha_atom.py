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

    def addLine(self, lineItem):
        pass

    def removeLine(self, lineItem):
        pass

    def itemChange(self, change, value):
        pass
