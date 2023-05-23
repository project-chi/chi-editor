import typing
from PyQt6 import QtGui
from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtWidgets import QGraphicsTextItem, QWidget, QLineEdit, QGraphicsItem, QGraphicsSimpleTextItem


class AnswerField(QGraphicsItem):

    rect: QRectF
    background_color: QtGui.QColor

    content: str
    answer_field_edit: QLineEdit

    def __init__(self, x, y, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.rect = QRectF(x, y, 100., 70.)
        self.background_color = QtGui.QColor("lightgray")

        self.content = "default"

        self.setAcceptHoverEvents(True)

    def boundingRect(self) -> QRectF:
        return self.rect

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        super().hoverEnterEvent(event)
        self.answer_field_edit = QLineEdit(self.content)
        self.scene().addWidget(self.answer_field_edit)
        self.answer_field_edit.setFocus()

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        super().hoverLeaveEvent(event)
        self.content = self.answer_field_edit.text()
        self.answer_field_edit.deleteLater()

    def paint(self, painter: QtGui.QPainter, option: 'QStyleOptionGraphicsItem', widget: typing.Optional[QWidget] = ...) -> None:
        painter.save()

        painter.setBrush(self.background_color)
        painter.drawRect(self.rect)

        painter.setBrush(QtGui.QColor("black"))
        painter.drawText(self.rect, Qt.AlignmentFlag.AlignCenter, self.content)

        painter.restore()
