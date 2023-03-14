from __future__ import annotations

from typing import TYPE_CHECKING

from PyQt6 import QtCore

if TYPE_CHECKING:
    from chi_editor.bases.molecule.molecule import Molecule

from PyQt6.QtGui import QPen, QBrush, QColor, QFont, QPainter
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QStyleOptionGraphicsItem,
    QGraphicsSceneMouseEvent,
    QGraphicsScene,
)
from PyQt6.QtCore import QRectF, Qt, QVariant

from .line import Line


class AlphaAtom(QGraphicsItem):
    background_pen: QPen = QPen(QColor("white"), 1)
    text_pen: QPen = QPen(QColor("black"), 10)
    text_font: QFont = QFont("Helvetica", 40)
    brush: QBrush = QBrush(QColor("white"))
    rect: QRectF = QRectF(0, 0, 50, 50)

    molecule: Molecule
    text: str
    lines: list[Line]

    def __init__(self, element: str, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        self.text = element
        self.lines = []
        self.setZValue(1)

        self.setFlag(self.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(self.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(self.GraphicsItemFlag.ItemSendsScenePositionChanges)

        from chi_editor.bases.molecule.molecule import Molecule

        self.molecule = Molecule(self)

    def get_adjacent_atoms(self) -> list:
        adjacent_atoms = []
        for line in self.lines:
            adjacent_atoms.append(
                line.vertex2 if line.vertex1 == self else line.vertex1
            )
        return adjacent_atoms

    def setPos(self, pos: QtCore.QPointF) -> None:
        super().setPos(pos)
        self.molecule.update_atoms()

    def remove(self):
        list_of_lines = list(self.lines)
        for line in list_of_lines:
            line.remove()

        self.molecule.remove_atom(self)

        self.scene().removeItem(self)

    def add_line(self, new_line: Line) -> bool:
        for existing in self.lines:
            if (existing.vertex1, existing.vertex2) == (
                new_line.vertex1,
                new_line.vertex2,
            ) or (existing.vertex2, existing.vertex1) == (
                new_line.vertex1,
                new_line.vertex2,
            ):
                # another line with the same control points already exists
                return False
        self.lines.append(new_line)
        return True

    def add_to_canvas(self, canvas: QGraphicsScene):
        canvas.addItem(self)
        self.molecule.anchor.add_to_canvas(canvas)

    def itemChange(
        self, change: QGraphicsItem.GraphicsItemChange, value: QVariant
    ) -> QVariant:
        for line in self.lines:
            line.update_pixmap(self)
        return super().itemChange(change, value)

    def boundingRect(self) -> QRectF:
        # adjust for boarder width
        adjust = self.background_pen.width() / 2
        return self.rect.adjusted(-adjust, -adjust, adjust, adjust)

    def paint(
        self, painter: QPainter, option: QStyleOptionGraphicsItem, widget=None
    ) -> None:
        # save + restore to reset pen and brush
        painter.save()
        painter.setPen(self.background_pen)
        painter.setBrush(self.brush)
        painter.drawEllipse(self.rect)
        painter.setPen(self.text_pen)
        painter.setFont(self.text_font)
        painter.drawText(self.rect, Qt.AlignmentFlag.AlignCenter, self.text)
        painter.restore()

    def mouseMoveEvent(self, event: "QGraphicsSceneMouseEvent") -> None:
        super().mouseMoveEvent(event)
        self.molecule.update_atoms()
