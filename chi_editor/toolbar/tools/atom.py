from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem

from ...bases.tool import Tool
from ...bases.alpha_atom import AlphaAtom


class Atom(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_atom = AlphaAtom('H')
        new_atom.setPos(event.scenePos())
        new_atom.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        new_atom.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)

        self.canvas.addItem(new_atom)

    @property
    def asset(self) -> str:
        return 'atom'

