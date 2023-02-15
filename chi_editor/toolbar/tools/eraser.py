from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ...bases.alpha_atom import AlphaAtom
from ...bases.line import Line
from ...bases.tool import Tool
from ...future.chemistry.atom import Atom


class Eraser(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            items = self.canvas.items(event.scenePos(), Qt.ItemSelectionMode.IntersectsItemShape)
            if not items:
                return super(Eraser, self).mouse_press_event(event)

            if isinstance(items[0], AlphaAtom):
                atom = items[0]
                for line in atom.lines:
                    line.vertex1.remove_line(line)
                    line.vertex2.remove_line(line)
                    self.canvas.removeItem(line)

            elif isinstance(items[0], Line):
                line = items[0]
                line.vertex1.remove_line(line)
                line.vertex2.remove_line(line)

            self.canvas.removeItem(items[0])
        else:
            self.canvas.clear()


    @property
    def asset(self) -> str:
        return 'eraser'
