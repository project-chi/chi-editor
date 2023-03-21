from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ...bases.tool import Tool
from ...canvas import Canvas


class Arrow(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        super(Canvas, self.canvas).mousePressEvent(event)

    def mouse_move_event(self, event: QGraphicsSceneMouseEvent) -> None:
        super(Canvas, self.canvas).mouseMoveEvent(event)

    def mouse_release_event(self, event: QGraphicsSceneMouseEvent) -> None:
        super(Canvas, self.canvas).mouseReleaseEvent(event)

    @property
    def asset(self) -> str:
        return "arrow"
