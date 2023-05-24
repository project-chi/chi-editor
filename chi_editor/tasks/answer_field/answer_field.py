from typing import Optional

from PyQt6 import QtGui
from PyQt6.QtCore import QRectF, Qt, QPointF
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import QWidget, QGraphicsItem, QGraphicsSceneHoverEvent, QStyleOptionGraphicsItem

from chi_editor.reactions.size_constants import Sizes

from chi_editor.tasks.answer_field.answer_field_menu import AnswerFieldMenu


class AnswerField(QGraphicsItem):
    rect: QRectF
    background_color: QtGui.QColor
    font: QFont

    content: str
    answer_field_menu: "AnswerFieldMenu"

    editable: bool

    def __init__(self, x, y, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.rect = QRectF(QPointF(x, y), Sizes.reagent_size)
        self.background_color = QtGui.QColor("lightgray")

        self.content = "default"

        self.setAcceptHoverEvents(True)
        self.answer_field_menu = AnswerFieldMenu(self, x, y)

        self.editable = True

    def boundingRect(self) -> QRectF:
        return self.rect

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        super().hoverEnterEvent(event)
        if self.editable:
            self.answer_field_menu.add_to_scene(self.scene())

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        super().hoverLeaveEvent(event)
        if self.editable:
            self.answer_field_menu.remove_from_scene(self.scene())

    def paint(self, painter: QtGui.QPainter, option: 'QStyleOptionGraphicsItem',
              widget: Optional[QWidget] = ...) -> None:
        painter.save()

        painter.setBrush(self.background_color)
        painter.drawRect(self.rect)

        painter.setBrush(QtGui.QColor("black"))
        painter.drawText(self.rect, Qt.AlignmentFlag.AlignCenter, self.content)

        painter.restore()

    def set_content(self, content: str) -> None:
        self.content = content
