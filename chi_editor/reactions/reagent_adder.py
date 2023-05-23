from typing import Optional, ClassVar

from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QGraphicsEllipseItem, QGraphicsSceneMouseEvent, QGraphicsItem, QStyleOptionGraphicsItem, \
    QWidget, QGraphicsSceneHoverEvent
from PyQt6.QtGui import QPainter, QColor


class ReagentAdder(QGraphicsEllipseItem):
    _reagent_list: list[QGraphicsItem]
    _item_highlighted: bool = False
    _highlight_color: QColor = QColor(0, 255, 50, 100)

    def __init__(self, reagent_list: list[QGraphicsItem], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._reagent_list = reagent_list
        self.setAcceptHoverEvents(True)

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent):
        self._reagent_list.append(QGraphicsItem())

    def paint(self, painter: QPainter, option: 'QStyleOptionGraphicsItem', widget: Optional[QWidget] = ...) -> None:
        painter.save()
        if self._item_highlighted:
            painter.setBrush(self._highlight_color)
            painter.drawEllipse(self.boundingRect())

        painter.drawLine(self._getTopMiddlePoint(), self._getBottomMiddlePoint())
        painter.drawLine(self._getRightMiddlePoint(), self._getLeftMiddlePoint())
        painter.restore()

    def _getTopMiddlePoint(self) -> QPointF:
        return self.boundingRect().topLeft() + QPointF(self.boundingRect().width() / 2, 0)

    def _getBottomMiddlePoint(self) -> QPointF:
        return self.boundingRect().bottomLeft() + QPointF(self.boundingRect().width() / 2, 0)

    def _getLeftMiddlePoint(self) -> QPointF:
        return self.boundingRect().topLeft() + QPointF(0, self.boundingRect().height() / 2)

    def _getRightMiddlePoint(self) -> QPointF:
        return self.boundingRect().topRight() + QPointF(0, self.boundingRect().height() / 2)

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        self._item_highlighted = True
        self.update()

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        self._item_highlighted = False
        self.update()
