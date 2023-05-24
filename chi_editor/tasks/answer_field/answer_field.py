import typing

from PyQt6 import QtGui
from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import QWidget, QGraphicsItem

from chi_editor.tasks.answer_field.answer_field_menu import AnswerFieldMenu


class AnswerField(QGraphicsItem):
    rect: QRectF
    background_color: QtGui.QColor
    font: QFont

    content: str
    answer_field_menu: AnswerFieldMenu

    editable: bool

    def __init__(self, x, y, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.rect = QRectF(x, y, 400, 300)
        self.background_color = QtGui.QColor("lightgray")

        self.content = "default"

        self.setAcceptHoverEvents(True)
        self.answer_field_menu = AnswerFieldMenu(self, x, y)

        self.font = QFont()
        self.font.setPointSizeF(40)

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
              widget: typing.Optional[QWidget] = ...) -> None:
        painter.save()

        painter.setBrush(self.background_color)
        painter.drawRect(self.rect)

        painter.setBrush(QtGui.QColor("black"))
        painter.setFont(self.font)
        painter.drawText(self.rect, Qt.AlignmentFlag.AlignCenter, self.content)

        painter.restore()
