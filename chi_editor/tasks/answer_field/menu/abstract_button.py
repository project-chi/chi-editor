import typing
from PyQt6.QtCore import QRectF
from PyQt6.QtGui import QColor, QPainter
from PyQt6.QtWidgets import QGraphicsItem, QWidget


class AbstractButton(QGraphicsItem):

    rect: QRectF
    answer_field: 'AnswerField'

    def __init__(self, answer_field: 'AnswerField', x: float, y: float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.answer_field = answer_field
        self.rect = QRectF(x, y, 50, 50)
        self.setFlag(self.GraphicsItemFlag.ItemIsSelectable)

    def boundingRect(self) -> QRectF:
        return self.rect


    def paint(self, painter: QPainter, option: 'QStyleOptionGraphicsItem', widget: typing.Optional[QWidget] = ...) -> None:
        painter.save()

        painter.setBrush(self.background_color)
        painter.drawRect(self.rect)

        painter.restore()
