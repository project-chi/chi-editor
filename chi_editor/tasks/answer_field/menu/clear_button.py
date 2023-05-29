from typing import TYPE_CHECKING, ClassVar

from PyQt6.QtWidgets import QGraphicsSceneMouseEvent
from PyQt6.QtGui import QColor, QImage

from chi_editor.constants import ASSETS
from chi_editor.tasks.answer_field.menu.abstract_button import GraphicsAbstractButton

if TYPE_CHECKING:
    from chi_editor.tasks.answer_field.answer_field import AnswerField


class ClearButtonGraphics(GraphicsAbstractButton):
    _answer_field: 'AnswerField'

    def __init__(self, answer_field: 'AnswerField', x: float, y: float, width: float, height: float, *args, **kwargs):
        super().__init__(x, y, width, height, *args, **kwargs)
        self._answer_field = answer_field
        self.picture = QImage(str(ASSETS / "cleaner.png"))

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        self._answer_field.set_content("")
