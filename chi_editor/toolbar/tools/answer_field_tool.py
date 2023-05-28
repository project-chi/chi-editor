from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from chi_editor.bases.tool import Tool

from chi_editor.tasks.answer_field.answer_field import AnswerField
from chi_editor.canvas import Canvas


class AnswerFieldTool(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        pos: QPointF = event.scenePos()
        answer_field: "AnswerField" = AnswerField(pos.x(), pos.y())
        self.canvas.addItem(answer_field)

    def mouse_move_event(self, event: QGraphicsSceneMouseEvent) -> None:
        super(Canvas, self.canvas).mouseMoveEvent(event)

    def mouse_release_event(self, event: QGraphicsSceneMouseEvent) -> None:
        super(Canvas, self.canvas).mouseReleaseEvent(event)

    @property
    def picture(self) -> str:
        return "arrow"
