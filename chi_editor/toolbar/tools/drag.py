from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsView

from ...bases.tool import Tool


class Drag(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        for view in self.canvas.views():
            view.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)

    def mouse_move_event(self, event: QGraphicsSceneMouseEvent) -> None:
        pass

    def mouse_release_event(self, event: QGraphicsSceneMouseEvent) -> None:
        for view in self.canvas.views():
            view.setDragMode(QGraphicsView.DragMode.NoDrag)

    @property
    def picture(self) -> str:
        return "hand"
