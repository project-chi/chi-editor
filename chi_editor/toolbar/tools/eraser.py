from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ...bases.alpha_atom import AlphaAtom
from ...bases.line import Line
from ...bases.tool import Tool


class Eraser(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            items = self.canvas.items(
                event.scenePos(), Qt.ItemSelectionMode.IntersectsItemShape
            )
            if not items:
                return super(Eraser, self).mouse_press_event(event)
            elif isinstance(items[0], AlphaAtom) or isinstance(items[0], Line):
                items[0].remove()
        else:
            self.canvas.clear()

    @property
    def picture(self) -> str:
        return "eraser"
