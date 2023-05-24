from typing import TYPE_CHECKING

from PyQt6.QtGui import QColor

from chi_editor.popup_canvas import PopupCanvas
from chi_editor.tasks.answer_field.menu.abstract_button import AbstractButton

if TYPE_CHECKING:
    from chi_editor.tasks.answer_field.answer_field import AnswerField


class EditButton(AbstractButton):
    background_color: QColor
    answer_window: PopupCanvas

    def __init__(self, answer_field: 'AnswerField', x: float, y: float, *args, **kwargs):
        super().__init__(answer_field, x, y, *args, **kwargs)
        self.background_color = QColor("gray")

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        self.answer_window = PopupCanvas()
        self.answer_window.setGeometry(100, 100, 400, 200)
        self.answer_window.canvas_changed.connect(self.get_input)
        self.answer_window.show()

    def get_input(self, answer: str):
        self.answer_field.set_content(answer)
        self.answer_window.close()
