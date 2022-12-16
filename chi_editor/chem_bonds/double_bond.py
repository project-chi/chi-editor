from PyQt6.QtGui import QPainter, QColor, QPen
from PyQt6.QtCore import QPointF

from ..bases.line import Line


class DoubleBond(Line):
    multiplicity: int = 2

    def paintLine(self, painter: QPainter) -> None:
        pen = QPen(QColor("black"), 3)
        painter.setPen(pen)

        # draw straight line
        painter.drawLine(QPointF(self.width / 3, 0),
                         QPointF(self.width / 3, self.height))
        painter.drawLine(QPointF(2 * self.width / 3, 0),
                         QPointF(2 * self.width / 3, self.height))
