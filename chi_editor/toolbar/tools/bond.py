from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsTextItem, QGraphicsItem

from ...bases.tool import Tool
from ...bases.line import Line


class Bond(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        selected = self.canvas.selectedItems()
        if len(selected) != 2:
            return

        line = Line(selected[0], selected[1])
        self.canvas.addItem(line)

    @property
    def asset(self) -> str:
        return 'bond'

