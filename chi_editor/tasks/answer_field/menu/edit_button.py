from typing import TYPE_CHECKING

from PyQt6.QtCore import QRectF
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QGraphicsItem

from chi_editor.tasks.answer_field.menu.abstract_button import AbstractButton

if TYPE_CHECKING:
    from chi_editor.tasks.answer_field.answer_field import AnswerField


class EditButton(AbstractButton):

    background_color: QColor

    def __init__(self, answer_field: 'AnswerField', x: float, y: float, *args, **kwargs):
        super().__init__(answer_field, x, y, *args, **kwargs)
        self.background_color = QColor("gray")

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        pass
