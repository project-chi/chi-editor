from typing import TYPE_CHECKING

from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QColor, QPen

from chi_editor.bases.line import Line

if TYPE_CHECKING:
    from PyQt6.QtGui import QPainter


class SingleBond(Line):
    multiplicity = 1

    def paint_line(self, painter: "QPainter") -> "None":
        pen = QPen(QColor("black"), 3)
        painter.setPen(pen)

        # draw straight line
        painter.drawLine(QPointF(self.width / 2, 0),
                         QPointF(self.width / 2, self.height))