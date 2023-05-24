from typing import TYPE_CHECKING, ClassVar

from PyQt6.QtCore import QRectF, QPointF, QSizeF
from PyQt6.QtWidgets import QGraphicsItem, QGraphicsScene

from chi_editor.tasks.answer_field.menu.clear_button import ClearButton
from chi_editor.tasks.answer_field.menu.edit_button import EditButton
from chi_editor.reactions.size_constants import Sizes

if TYPE_CHECKING:
    from chi_editor.tasks.answer_field.answer_field import AnswerField


class AnswerFieldMenu:
    _reagent_menu_size: ClassVar[QSizeF] = QSizeF(Sizes.reagent_size.width() * 0.1, Sizes.reagent_size.height() * 0.2)
    _reagent_menu_button_size: ClassVar[QSizeF] = QSizeF(_reagent_menu_size.width(), _reagent_menu_size.height() / 2)
    answer_field: 'AnswerField'
    rect: QRectF

    clear_button: QGraphicsItem
    edit_button: QGraphicsItem

    def __init__(self, answer_field: 'AnswerField', x: float, y: float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.rect = QRectF(QPointF(x, y), self._reagent_menu_size)
        self.answer_field = answer_field

        self.clear_button = ClearButton(answer_field, x, y, self._reagent_menu_button_size.width(),
                                        self._reagent_menu_button_size.height())
        self.edit_button = EditButton(answer_field, x, y + self._reagent_menu_button_size.height(),
                                      self._reagent_menu_button_size.width(), self._reagent_menu_button_size.height())

    def add_to_scene(self, scene: QGraphicsScene):
        scene.addItem(self.clear_button)
        scene.addItem(self.edit_button)

    def remove_from_scene(self, scene: QGraphicsScene):
        scene.removeItem(self.clear_button)
        scene.removeItem(self.edit_button)
