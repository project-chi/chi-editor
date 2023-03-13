from PyQt6.QtCore import QRectF, QPointF
from PyQt6.QtGui import QPen, QColor, QBrush, QPainter
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QStyleOptionGraphicsItem,
    QGraphicsSceneMouseEvent,
)

from chi_editor.bases.alpha_atom import AlphaAtom


def get_geometrical_center(atoms: list[AlphaAtom]) -> QPointF:
    position: QPointF = QPointF(0, 0)
    for atom in atoms:
        position = position + atom.pos()
    atoms_count: int = len(atoms)
    return position / atoms_count


class MoleculeDrawer(QGraphicsItem):
    background_pen: QPen = QPen(QColor("red"), 1)
    brush: QBrush = QBrush(QColor("red"))
    rect: QRectF = QRectF(0, 0, 20, 20)

    atoms: list[AlphaAtom]

    def __init__(self, atoms, *args, **kwargs):
        super().__init__()
        self.atoms = atoms
        self.setZValue(2)
        self.setFlag(self.GraphicsItemFlag.ItemSendsScenePositionChanges)
        self.setFlag(self.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(self.GraphicsItemFlag.ItemIsMovable)

    def update_position(self, atoms: list[AlphaAtom]) -> None:
        self.setPos(get_geometrical_center(atoms))

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
        painter.drawRect(self.rect)
        painter.restore()

    def mousePressEvent(self, event: "QGraphicsSceneMouseEvent") -> None:
        super().mousePressEvent(event)
        self.dpos: dict = {}
        for atom in self.atoms:
            self.dpos[atom] = atom.pos() - self.pos()

    def mouseMoveEvent(self, event: "QGraphicsSceneMouseEvent") -> None:
        super().mouseMoveEvent(event)
        for atom in self.atoms:
            atom.setPos(self.pos() + self.dpos[atom])
