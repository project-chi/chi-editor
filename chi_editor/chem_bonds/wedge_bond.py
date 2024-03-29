from typing import TYPE_CHECKING

from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QBrush, QColor, QPainterPath, QPen

from chi_editor.bases.line import Line

if TYPE_CHECKING:
    from PyQt6.QtGui import QPainter


class WedgeBond(Line):
    multiplicity = 1

    def paint_line(self, painter: "QPainter") -> "None":
        pen = QPen(QColor("black"), 0)
        brush = QBrush(QColor("black"))
        painter.setPen(pen)
        painter.setBrush(brush)

        wedge_path = QPainterPath(QPointF(0, 0))
        wedge_path.lineTo(self.width / 2, self.height)
        wedge_path.lineTo(self.width, 0)

        # draw straight line
        painter.drawPath(wedge_path)
