from PyQt6.QtGui import QPainter, QColor, QPen, QBrush, QPainterPath
from PyQt6.QtCore import QPointF

from ..bases.line import Line


class WedgeBond(Line):
    multiplicity = 1

    def paint_line(self, painter: QPainter) -> None:
        pen = QPen(QColor("black"), 0)
        brush = QBrush(QColor("black"))
        painter.setPen(pen)
        painter.setBrush(brush)

        wedge_path = QPainterPath(QPointF(0, 0))
        wedge_path.lineTo(self.width / 2, self.height)
        wedge_path.lineTo(self.width, 0)

        # draw straight line
        painter.drawPath(wedge_path)
