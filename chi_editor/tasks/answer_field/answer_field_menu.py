from typing import TYPE_CHECKING, ClassVar

from PyQt6.QtCore import QRectF, QPointF, QSizeF
from PyQt6.QtWidgets import QGraphicsItem, QGraphicsSceneMouseEvent, QGraphicsItemGroup

from chi_editor.tasks.answer_field.menu.clear_button import ClearButtonGraphics
from chi_editor.tasks.answer_field.menu.edit_button import EditButtonGraphics
from chi_editor.reactions.size_constants import Sizes

if TYPE_CHECKING:
    from chi_editor.tasks.answer_field.answer_field import AnswerField


class AnswerFieldMenu(QGraphicsItemGroup):
    _reagent_menu_size: ClassVar[QSizeF] = QSizeF(Sizes.reagent_size.width() * 0.4, Sizes.reagent_size.height() * 0.8)
    _reagent_menu_button_size: ClassVar[QSizeF] = QSizeF(_reagent_menu_size.width(), _reagent_menu_size.height() / 2)
    answer_field: 'AnswerField'

    clear_button: QGraphicsItem
    edit_button: QGraphicsItem

    def __init__(self, answer_field: 'AnswerField', *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.answer_field = answer_field
        boarder_offset = self.answer_field.pen.widthF()
        self.setPos(boarder_offset, boarder_offset)

        # setting self as parent is the same as adding to the group
        self.clear_button = ClearButtonGraphics(answer_field, 0, 0, self._reagent_menu_button_size.width(),
                                                self._reagent_menu_button_size.height(), parent=self)
        self.clear_button.hide()

        self.edit_button = EditButtonGraphics(answer_field, 0, 0 + self._reagent_menu_button_size.height(),
                                              self._reagent_menu_button_size.width(),
                                              self._reagent_menu_button_size.height(),
                                              parent=self)
        self.edit_button.hide()

    def show(self):
        self.clear_button.show()
        self.edit_button.show()

    def hide(self):
        self.clear_button.hide()
        self.edit_button.hide()

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        if self.clear_button.sceneBoundingRect().contains(event.scenePos()):
            self.clear_button.mousePressEvent(event)
        if self.edit_button.sceneBoundingRect().contains(event.scenePos()):
            self.edit_button.mousePressEvent(event)
