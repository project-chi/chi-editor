from typing import Optional, ClassVar

from PyQt6 import QtGui
from PyQt6.QtCore import QRectF, Qt, QPointF
from PyQt6.QtGui import QFont, QColor, QPen
from PyQt6.QtWidgets import QWidget, QGraphicsItem, QGraphicsSceneHoverEvent, QStyleOptionGraphicsItem, \
    QGraphicsSceneMouseEvent

from chi_editor.tasks.tasks_size_constants import Sizes

from chi_editor.tasks.answer_field.answer_field_menu import AnswerFieldMenu


class AnswerField(QGraphicsItem):
    background_color: QColor
    font: QFont
    pen: ClassVar[QPen] = QPen(QColor("black"), Sizes.reagent_boarder_width)

    content: str
    answer_field_menu: "AnswerFieldMenu"

    editable: bool

    def __init__(self, x: float, y: float, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setPos(x, y)

        self.background_color = QtGui.QColor("lightgray")
        self.content = "c"
        self.editable = True

        self.answer_field_menu = AnswerFieldMenu(self, parent=self)

        self.setAcceptHoverEvents(True)
        self.setFiltersChildEvents(True)

    def boundingRect(self) -> QRectF:
        boarder_width = self.pen.widthF()
        bounding_rect = QRectF(QPointF(0, 0), Sizes.reagent_size)
        bounding_rect.adjust(
            boarder_width, boarder_width,
            boarder_width, boarder_width
        )
        return bounding_rect

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        if self.editable:
            self.answer_field_menu.show()

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        if self.editable:
            self.answer_field_menu.hide()

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        self.answer_field_menu.mousePressEvent(event)

    def paint(self, painter: QtGui.QPainter, option: 'QStyleOptionGraphicsItem',
              widget: Optional[QWidget] = ...) -> None:
        painter.save()

        painter.setPen(self.pen)
        painter.setBrush(self.background_color)
        painter.drawRect(self.boundingRect())

        painter.setBrush(QtGui.QColor("black"))
        painter.drawText(self.boundingRect(), Qt.AlignmentFlag.AlignCenter, self.content)

        painter.restore()

    def set_content(self, content: str) -> None:
        self.content = content
        self.update()
