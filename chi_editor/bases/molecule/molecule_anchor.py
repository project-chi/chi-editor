import weakref

from PyQt6.QtCore import QRectF, QPointF
from PyQt6.QtGui import QPen, QColor, QBrush, QPainter
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QStyleOptionGraphicsItem,
    QGraphicsSceneMouseEvent,
    QGraphicsScene,
)

from chi_editor.bases.alpha_atom import AlphaAtom


def get_geometrical_center(atoms: list[AlphaAtom]) -> QPointF:
    position: QPointF = QPointF(0, 0)
    for atom in atoms:
        position = position + atom.pos()
    atoms_count: int = len(atoms)
    return position / atoms_count


class MoleculeAnchor(QGraphicsItem):
    background_pen: QPen = QPen(QColor("red"), 1)
    brush: QBrush = QBrush(QColor("red"))
    rect: QRectF = QRectF(-10, -10, 20, 20)

    atoms: weakref.WeakSet[AlphaAtom]
    atom_delta_positions: dict

    def __init__(self, atoms: weakref.WeakSet[AlphaAtom], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.atoms = atoms
        self.setZValue(2)
        self.setFlag(self.GraphicsItemFlag.ItemSendsScenePositionChanges)
        self.setFlag(self.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(self.GraphicsItemFlag.ItemIsMovable)

    def update_position(self, atoms: list[AlphaAtom]) -> None:
        self.setPos(get_geometrical_center(atoms) + AlphaAtom.rect.center())

    def boundingRect(self) -> QRectF:
        # adjust for boarder width
        adjust = self.background_pen.width() / 2
        return self.rect.adjusted(-adjust, -adjust, adjust, adjust)

    def add_to_canvas(self, canvas: QGraphicsScene):
        canvas.addItem(self)

    def remove(self):
        if self.scene() is not None:
            self.scene().removeItem(self)

    def paint(
        self, painter: QPainter, option: QStyleOptionGraphicsItem, widget=None
    ) -> None:
        # save + restore to reset pen and brush
        painter.save()
        painter.setPen(self.background_pen)
        painter.setBrush(self.brush)
        painter.drawRect(self.rect)
        painter.restore()

    def mouseMoveEvent(self, event: "QGraphicsSceneMouseEvent") -> None:
        super().mouseMoveEvent(event)
        for atom in self.atoms:
            atom.moveBy(
                event.pos().x() - event.lastPos().x(),
                event.pos().y() - event.lastPos().y(),
            )
