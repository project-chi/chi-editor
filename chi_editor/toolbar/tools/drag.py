from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ...bases.tool import Tool


class Drag(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        pass

    def mouse_move_event(self, event: QGraphicsSceneMouseEvent) -> None:
        pass

    def mouse_release_event(self, event: QGraphicsSceneMouseEvent) -> None:
        pass

    @property
    def asset(self) -> str:
        return 'arrow'
