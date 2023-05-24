from typing import TYPE_CHECKING

from PyQt6.QtGui import QColor

from chi_editor.tasks.answer_field.menu.abstract_button import AbstractButton

if TYPE_CHECKING:
    from chi_editor.tasks.answer_field.answer_field import AnswerField


class ClearButton(AbstractButton):
    background_color: QColor

    def __init__(self, answer_field: 'AnswerField', x: float, y: float, width: float, height: float, *args, **kwargs):
        super().__init__(answer_field, x, y, width, height, *args, **kwargs)
        self.background_color = QColor("red")

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        super().mousePressEvent(event)
        self.answer_field.content = ""
