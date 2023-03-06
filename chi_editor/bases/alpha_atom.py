from typing import List

from PyQt6.QtGui import QPen, QBrush, QColor, QFont, QPainter
from PyQt6.QtWidgets import QGraphicsItem, QStyleOptionGraphicsItem
from PyQt6.QtCore import QRectF, Qt, QVariant

from .line import Line


class AlphaAtom(QGraphicsItem):
    background_pen: QPen = QPen(QColor("white"), 1)
    text_pen: QPen = QPen(QColor("black"), 10)
    text_font: QFont = QFont("Helvetica", 40)
    brush: QBrush = QBrush(QColor("white"))
    rect: QRectF = QRectF(0, 0, 50, 50)

    _text: str
    _lines: list[Line]

    def __init__(self, element: str, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.text = element
        self.lines = []
        self.setZValue(1)
        self.setFlags(self.GraphicsItemFlag.ItemSendsScenePositionChanges)

    def get_adjacent_atoms(self) -> list:
        adjacent_atoms = []
        for line in self.lines:
            adjacent_atoms.append(line.vertex2 if line.vertex1 == self else line.vertex1)
        return adjacent_atoms

    def get_molecule_atoms(self) -> list[QGraphicsItem]:
        marked: list[AlphaAtom] = []
        queue: list[AlphaAtom] = [self]
        while queue:
            current_atom: AlphaAtom = queue.pop()
            marked.append(current_atom)
            queue.extend([x for x in current_atom.get_adjacent_atoms() if x not in marked and x not in queue])
        return marked

    def remove(self):
        list_of_lines = list(self.lines)
        for line in list_of_lines:
            line.remove()
        self.scene().removeItem(self)

    def add_line(self, new_line: Line) -> bool:
        for existing in self.lines:
            if (existing.vertex1, existing.vertex2) == (new_line.vertex1, new_line.vertex2) \
                    or (existing.vertex2, existing.vertex1) == (new_line.vertex1, new_line.vertex2):
                # another line with the same control points already exists
                return False
        self.lines.append(new_line)
        return True

    def itemChange(self, change: QGraphicsItem.GraphicsItemChange, value: QVariant) -> QVariant:
        for line in self.lines:
            line.update_pixmap(self)
        return super().itemChange(change, value)

    def boundingRect(self) -> QRectF:
        # adjust for boarder width
        adjust = self.background_pen.width() / 2
        return self.rect.adjusted(-adjust, -adjust, adjust, adjust)

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget=None) -> None:
        # save + restore to reset pen and brush
        painter.save()
        painter.setPen(self.background_pen)
        painter.setBrush(self.brush)
        painter.drawEllipse(self.rect)
        painter.setPen(self.text_pen)
        painter.setFont(self.text_font)
        painter.drawText(self.rect, Qt.AlignmentFlag.AlignCenter, self.text)
        painter.restore()
