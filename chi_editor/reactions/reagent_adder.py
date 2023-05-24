from typing import Optional, ClassVar
from enum import Enum

from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QGraphicsEllipseItem, QGraphicsSceneMouseEvent, QGraphicsItem, QStyleOptionGraphicsItem, \
    QWidget, QGraphicsSceneHoverEvent, QGraphicsRectItem
from PyQt6.QtGui import QPainter, QColor

from chi_editor.reactions.size_constants import Sizes


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
        next_pos = self._getNextItemPos()
        new_item = QGraphicsRectItem(next_pos.x(), next_pos.y(), Sizes.reagent_size.width(),
                                     Sizes.reagent_size.height())

        self._reagent_list.append(new_item)
        self.setPos(self._getNextAddPos())

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

    def _getNextAddPos(self) -> QPointF:
        match self._growth_direction:
            case GrowthDirection.LEFT:
                return self.pos() + QPointF(
                    -1 * (Sizes.default_gap + Sizes.plus_size.width() + Sizes.default_gap + Sizes.reagent_size.width()),
                    0)
            case GrowthDirection.RIGHT:
                return self.pos() + QPointF(
                    Sizes.reagent_size.width() + Sizes.default_gap + Sizes.plus_size.width() + Sizes.default_gap, 0)

    def _getNextItemPos(self) -> QPointF:
        match self._growth_direction:
            case GrowthDirection.LEFT:
                return self._getLastItemPos() + QPointF(
                    -1 * (Sizes.default_gap + Sizes.plus_size.width() + Sizes.default_gap + Sizes.reagent_size.width()),
                    0)
            case GrowthDirection.RIGHT:
                return self._getLastItemPos() + QPointF(
                    Sizes.reagent_size.width() + Sizes.default_gap + Sizes.plus_size.width() + Sizes.default_gap, 0)

    def _getLastItemPos(self) -> QPointF:
        item = self._reagent_list.pop()
        self._reagent_list.append(item)
        return item.sceneBoundingRect().topLeft()

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        self._item_highlighted = True
        self.update()

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        self._item_highlighted = False
        self.update()
