from typing import TYPE_CHECKING

from PyQt6.QtWidgets import QGraphicsSceneMouseEvent
from PyQt6.QtGui import QColor

from chi_editor.popup_canvas import PopupCanvas
from chi_editor.tasks.answer_field.menu.abstract_button import GraphicsAbstractButton

if TYPE_CHECKING:
    from chi_editor.tasks.answer_field.answer_field import AnswerField


class EditButtonGraphics(GraphicsAbstractButton):
    background_color: QColor
    answer_window: PopupCanvas
    _answer_field: 'AnswerField'

    def __init__(self, answer_field: 'AnswerField', x: float, y: float, width: float, height: float, *args, **kwargs):
        super().__init__(x, y, width, height, *args, **kwargs)
        self._answer_field = answer_field
        self.background_color = QColor("gray")

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        self.answer_window = PopupCanvas()
        self.answer_window.setGeometry(100, 100, 400, 200)
        self.answer_window.canvas_changed.connect(self.get_input)
        self.answer_window.show()

    def get_input(self, answer: str):
        self._answer_field.set_content(answer)
        self.answer_window.close()
