from PyQt6.QtGui import QPen, QBrush, QColor, QFont
from PyQt6.QtWidgets import QGraphicsItem
from PyQt6.QtCore import QRectF, Qt, QVariant

from .line import Line


class AlphaAtom(QGraphicsItem):
    text: str
    pen: QPen = QPen(QColor("black"), 0)
    brush: QBrush = QBrush(QColor("black"))
    rect: QRectF = QRectF(0, 0, 100, 100)
    lines: list[Line] = []

    def __init__(self, element: str, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.text = element
        self.setFlags(self.GraphicsItemFlag.ItemSendsScenePositionChanges)

    def addLine(self, newLine: Line) -> bool:
        for existing in self.lines:
            if (existing.vertex1, existing.vertex2) == (newLine.vertex1, newLine.vertex2):
                # another line with the same control points already exists
                return False
        self.lines.append(newLine)
        return True

    def removeLine(self, line: Line) -> bool:
        for existing in self.lines:
            if (existing.vertex1, existing.vertex2) == (line.vertex1, line.vertex2):
                self.scene().removeItem(existing)
                self.lines.remove(existing)
                return True
        return False

    def itemChange(self, change: QGraphicsItem.GraphicsItemChange, value: QVariant) -> QVariant:
        for line in self.lines:
            line.update_pixmap(self)
        return super().itemChange(change, value)

    def boundingRect(self) -> QRectF:
        # adjust for boarder width
        adjust = self.pen.width() / 2
        return self.rect.adjusted(-adjust, -adjust, adjust, adjust)

    def paint(self, painter, option, widget=None) -> None:
        # save + restore to reset pen and brush
        painter.save()
        painter.setPen(self.pen)
        painter.setBrush(self.brush)
        painter.drawEllipse(self.rect)
        text_pen = QPen(QColor("white"), 10)
        painter.setPen(text_pen)
        painter.setFont(QFont("Impact", 30))
        painter.drawText(self.rect, Qt.AlignmentFlag.AlignCenter, self.text)
        painter.restore()
