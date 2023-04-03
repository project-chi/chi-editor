from typing import TYPE_CHECKING
from weakref import WeakSet

from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtGui import QBrush, QColor, QImage, QPainter, QPen
from PyQt6.QtWidgets import QGraphicsItem

from chi_editor.bases.sources import BASIC_RECTANGLE
from chi_editor.constants import ASSETS

if TYPE_CHECKING:
    from typing import ClassVar

    from PyQt6.QtWidgets import QGraphicsScene, QGraphicsSceneMouseEvent

    from chi_editor.bases.alpha_atom import AlphaAtom


def get_geometrical_center(atoms: "list[AlphaAtom]") -> "QPointF":
    position: "QPointF" = QPointF(0, 0)
    for atom in atoms:
        position = position + atom.pos()
    atoms_count: int = len(atoms)
    return position / atoms_count


class MoleculeAnchor(QGraphicsItem):
    picture: "ClassVar[QImage]" = QImage(str(ASSETS / "anchor.png"))

    background_pen: "ClassVar[QPen]" = QPen(QColor("lightgray"), 1)
    brush: "ClassVar[QBrush]" = QBrush(QColor("lightgray"))
    rect: "ClassVar[QRectF]" = QRectF(
        BASIC_RECTANGLE.center().x() - 8, BASIC_RECTANGLE.center().y() - 8, 16, 16
    )

    atoms: "WeakSet[AlphaAtom]"

    def __init__(self, atoms: "WeakSet[AlphaAtom]", *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.atoms = atoms
        self.setZValue(2)
        self.setFlag(self.GraphicsItemFlag.ItemSendsScenePositionChanges)
        self.setFlag(self.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(self.GraphicsItemFlag.ItemIsMovable)

    def update_position(self, atoms: "list[AlphaAtom]") -> "None":
        self.setPos(get_geometrical_center(atoms))

    def boundingRect(self) -> "QRectF":
        # adjust for boarder width
        adjust = self.background_pen.width() / 2
        return self.rect.adjusted(-adjust, -adjust, adjust, adjust)

    def add_to_canvas(self, canvas: "QGraphicsScene") -> "None":
        canvas.addItem(self)

    def remove(self) -> "None":
        if self.scene() is not None:
            self.scene().removeItem(self)

    def paint(self, painter: "QPainter", *_) -> "None":
        # save + restore to reset pen and brush
        painter.save()
        painter.setPen(self.background_pen)
        painter.setBrush(self.brush)
        painter.drawEllipse(self.rect)
        painter.drawImage(self.rect, self.picture)
        painter.restore()

    def mouseMoveEvent(self, event: "QGraphicsSceneMouseEvent") -> "None":
        super().mouseMoveEvent(event)
        for atom in self.atoms:
            atom.moveBy(
                event.pos().x() - event.lastPos().x(),
                event.pos().y() - event.lastPos().y(),
            )
