from typing import TYPE_CHECKING

from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QColor, QPen

from chi_editor.bases.line import Line

if TYPE_CHECKING:
    from PyQt6.QtGui import QPainter


class DoubleBond(Line):
    multiplicity = 2

    def paint_line(self, painter: "QPainter") -> "None":
        pen = QPen(QColor("black"), 3)
        painter.setPen(pen)

        # draw straight line
        painter.drawLine(
            QPointF(self.width / 3, 0),
            QPointF(self.width / 3, self.height)
        )
        painter.drawLine(
            QPointF(2 * self.width / 3, 0),
            QPointF(2 * self.width / 3, self.height)
        )
