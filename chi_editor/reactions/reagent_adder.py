from typing import Optional, ClassVar
from enum import Enum

from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QGraphicsEllipseItem, QGraphicsSceneMouseEvent, QGraphicsItem, QStyleOptionGraphicsItem, \
    QWidget, QGraphicsSceneHoverEvent, QGraphicsRectItem
from PyQt6.QtGui import QPainter, QColor


class GrowthDirection(Enum):
    LEFT = 1
    RIGHT = 2


class ReagentAdder(QGraphicsEllipseItem):
    _reagent_list: list[QGraphicsItem]
    _item_highlighted: bool = False
    _highlight_color: ClassVar[QColor] = QColor(0, 255, 50, 100)
    _growth_direction: GrowthDirection

    def __init__(self, reagent_list: list[QGraphicsItem], growth_direction: GrowthDirection, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._reagent_list = reagent_list
        self._growth_direction = growth_direction
        self.setAcceptHoverEvents(True)

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent):
        new_item = QGraphicsRectItem(len(self._reagent_list))

        self._reagent_list.append()

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

    # def _getNextPos(self) -> QPointF:
    #     match self._growth_direction:
    #         case Direction.LEFT:
    #             return QPointF(-1 * len(self._reagent_list) * ())

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        self._item_highlighted = True
        self.update()

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        self._item_highlighted = False
        self.update()
