import typing
from typing import ClassVar

from PyQt6.QtCore import QRectF
from PyQt6.QtGui import QPainter, QColor, QImage
from PyQt6.QtWidgets import QGraphicsItem, QWidget, QStyleOptionGraphicsItem


class GraphicsAbstractButton(QGraphicsItem):
    _bounding_rect: QRectF
    picture: "ClassVar[QImage]"

    def __init__(self, x: float, y: float, width: float, height: float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._bounding_rect = QRectF(x, y, width, height)
        self.setFlag(self.GraphicsItemFlag.ItemIsSelectable)

    def boundingRect(self) -> QRectF:
        return self._bounding_rect

    def paint(self, painter: QPainter, option: 'QStyleOptionGraphicsItem',
              widget: typing.Optional[QWidget] = ...) -> None:
        painter.save()

        painter.drawImage(self._bounding_rect, self.picture)

        painter.restore()
