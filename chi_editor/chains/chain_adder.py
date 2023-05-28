from typing import Optional, ClassVar, TYPE_CHECKING
from enum import Enum

from PyQt6.QtCore import QPointF, QObject, QEvent, QRectF
from PyQt6.QtWidgets import QGraphicsEllipseItem, QGraphicsSceneMouseEvent, QGraphicsItem, QStyleOptionGraphicsItem, \
    QWidget, QGraphicsSceneHoverEvent
from PyQt6.QtGui import QPainter, QColor

from chi_editor.tasks.tasks_size_constants import Sizes

from chi_editor.tasks.answer_field.answer_field import AnswerField

if TYPE_CHECKING:
    from chi_editor.chains.chain_item_group import ChainItemGroup


class GrowthDirection(Enum):
    LEFT = 1
    RIGHT = 2


class ChainReagentAdder(QGraphicsEllipseItem):
    _chain_group: 'ChainItemGroup'
    _reagent_list: list[QGraphicsItem]
    _item_highlighted: bool = False
    _highlight_color: ClassVar[QColor] = QColor(0, 255, 50, 100)
    _growth_direction: GrowthDirection

    def __init__(self, reaction_group: 'ChainItemGroup', reagent_list: list[QGraphicsItem],
                 growth_direction: GrowthDirection, pos: QPointF, *args,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.setPos(pos)
        self.setRect(QRectF(QPointF(0, 0), Sizes.add_item_size))
        self._chain_group = reaction_group
        self._reagent_list = reagent_list
        self._growth_direction = growth_direction
        self.setAcceptHoverEvents(True)

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent):
        next_pos = self._getNextItemPos()
        new_item = AnswerField(next_pos.x(), next_pos.y())

        # add new item in corresponding position
        match self._growth_direction:
            case GrowthDirection.LEFT:
                self._reagent_list.insert(0, new_item)
            case GrowthDirection.RIGHT:
                self._reagent_list.append(new_item)

        self._chain_group.addToGroup(new_item)
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
                    -1 * (
                                Sizes.default_gap + Sizes.arrow_size.width() + Sizes.default_gap + Sizes.reagent_size.width()),
                    0)
            case GrowthDirection.RIGHT:
                return self.pos() + QPointF(
                    Sizes.reagent_size.width() + Sizes.default_gap + Sizes.arrow_size.width() + Sizes.default_gap, 0)

    def _getNextItemPos(self) -> QPointF:
        match self._growth_direction:
            case GrowthDirection.LEFT:
                return self._getFirstItemPos() + QPointF(
                    -1 * (
                                Sizes.default_gap + Sizes.arrow_size.width() + Sizes.default_gap + Sizes.reagent_size.width()),
                    0) - QPointF(Sizes.reagent_boarder_width, Sizes.reagent_boarder_width)
            case GrowthDirection.RIGHT:
                return self._getLastItemPos() + QPointF(
                    Sizes.reagent_size.width() + Sizes.default_gap + Sizes.arrow_size.width() + Sizes.default_gap,
                    0) - QPointF(Sizes.reagent_boarder_width, Sizes.reagent_boarder_width)

    def _getLastItemPos(self) -> QPointF:
        return self._reagent_list[-1].sceneBoundingRect().topLeft()

    def _getFirstItemPos(self) -> QPointF:
        return self._reagent_list[0].sceneBoundingRect().topLeft()

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        self._item_highlighted = True
        self.update()

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        self._item_highlighted = False
        self.update()
