from typing import TYPE_CHECKING

from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QColor, QPen

from chi_editor.bases.line import Line

if TYPE_CHECKING:
    from PyQt6.QtGui import QPainter


class TripleBond(Line):
    multiplicity = 3

    def paint_line(self, painter: "QPainter") -> "None":
        pen = QPen(QColor("black"), 3)
        painter.setPen(pen)

        # draw straight line
        painter.drawLine(QPointF(self.width / 4, 0),
                         QPointF(self.width / 4, self.height))
        painter.drawLine(QPointF(self.width / 2, 0),
                         QPointF(self.width / 2, self.height))
        painter.drawLine(QPointF(3 * self.width / 4, 0),
                         QPointF(3 * self.width / 4, self.height))
