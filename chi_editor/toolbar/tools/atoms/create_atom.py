from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from chi_editor.bases.alpha_atom import AlphaAtom
from chi_editor.bases.tool import Tool


class CreateAtom(Tool):
    _element: str

    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_atom = AlphaAtom(self._element)
        new_atom.setPos(event.scenePos() - new_atom.sceneBoundingRect().center())
        new_atom.add_to_canvas(self.canvas)

    @property
    def picture(self) -> str:
        return "carbon"
