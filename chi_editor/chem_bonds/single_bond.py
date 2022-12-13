from PyQt6.QtGui import QPainter, QColor, QPen
from PyQt6.QtCore import QPointF

from ..bases.line import Line


class SingleBond(Line):
    multiplicity = 1

    def paintLine(self, painter: QPainter) -> None:
        pen = QPen(QColor("black"), 3)
        painter.setPen(pen)

        # draw straight line
        painter.drawLine(QPointF(self.width / 2, 0),
                         QPointF(self.width / 2, self.height))