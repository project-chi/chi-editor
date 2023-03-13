from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ...bases.alpha_atom import AlphaAtom
from ...bases.tool import Tool


class Atom(Tool):
    _element: str

    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_atom = AlphaAtom(self._element)
        new_atom.setPos(event.scenePos() - new_atom.sceneBoundingRect().center())
        new_atom.add_to_canvas(self.canvas)

    @property
    def asset(self) -> str:
        return "carbon"
