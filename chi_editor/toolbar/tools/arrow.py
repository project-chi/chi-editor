from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ...bases.tool import Tool
from ...canvas import Canvas


class Arrow(Tool):
    def action(self, event: QGraphicsSceneMouseEvent) -> None:
        super(Canvas, self.canvas).mousePressEvent(event)

    @property
    def asset(self) -> str:
        return 'arrow'
