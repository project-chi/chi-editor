from typing import TYPE_CHECKING

from PyQt6.QtCore import QRectF
from PyQt6.QtWidgets import QGraphicsItem, QGraphicsScene

from chi_editor.tasks.answer_field.menu.clear_button import ClearButton
from chi_editor.tasks.answer_field.menu.edit_button import EditButton

if TYPE_CHECKING:
    from chi_editor.tasks.answer_field.answer_field import AnswerField


class AnswerFieldMenu:
    answer_field: 'AnswerField'
    rect: QRectF

    clear_button: QGraphicsItem
    edit_button: QGraphicsItem

    def __init__(self, answer_field: 'AnswerField', x: float, y: float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.rect = QRectF(x, y, 50, 100)
        self.answer_field = answer_field

        self.clear_button = ClearButton(answer_field, x, y)
        self.edit_button = EditButton(answer_field, x, y + 50)

    def add_to_scene(self, scene: QGraphicsScene):
        scene.addItem(self.clear_button)
        scene.addItem(self.edit_button)

    def remove_from_scene(self, scene: QGraphicsScene):
        scene.removeItem(self.clear_button)
        scene.removeItem(self.edit_button)
