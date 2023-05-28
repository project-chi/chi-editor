from typing import Optional

from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtGui import QPainter, QPainterPath
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QGraphicsItemGroup,
    QGraphicsEllipseItem,
    QGraphicsSceneHoverEvent,
    QStyleOptionGraphicsItem,
    QWidget,
    QGraphicsSceneMouseEvent
)

from chi_editor.chains.chain_adder import ChainReagentAdder, GrowthDirection
from chi_editor.tasks.tasks_size_constants import Sizes
from chi_editor.tasks.answer_field.answer_field import AnswerField


class ChainItemGroup(QGraphicsItemGroup):
    _reagent_items: list[QGraphicsItem]
    _add_left_item: QGraphicsEllipseItem
    _add_right_item: QGraphicsEllipseItem

    def __init__(self, x: float, y: float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setPos(x, y)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFiltersChildEvents(False)
        self.setAcceptHoverEvents(False)

        self._reagent_items = []
        self._addAddButtons()
        self._addInitialItems()

    def boundingRect(self) -> QRectF:
        return self.childrenBoundingRect()

    def _addAddButtons(self) -> None:
        add_left_top_left_point = self.pos() + QPointF(-1 * (
                Sizes.arrow_size.width() / 2 + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap + Sizes.add_item_size.width()),
                                                       -1 * Sizes.add_item_size.height() / 2)
        self._add_left_item = ChainReagentAdder(self, self._reagent_items, GrowthDirection.LEFT,
                                                add_left_top_left_point)

        add_right_top_left_point = self.pos() + QPointF(
            Sizes.arrow_size.width() / 2 + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap,
            -1 * Sizes.add_item_size.height() / 2)
        self._add_right_item = ChainReagentAdder(self, self._reagent_items, GrowthDirection.RIGHT,
                                                 add_right_top_left_point)

        self.addToGroup(self._add_right_item)
        self.addToGroup(self._add_left_item)

    def _addInitialItems(self) -> None:
        first_reagent_point = self.pos() + QPointF(-1 * (Sizes.side_items_offset + Sizes.reagent_size.width()),
                                                   -1 * Sizes.reagent_size.height() / 2)
        first_reagent = AnswerField(first_reagent_point.x(), first_reagent_point.y())
        self._reagent_items.append(first_reagent)
        self.addToGroup(first_reagent)

        second_reagent_point = self.pos() + QPointF(Sizes.side_items_offset, -1 * Sizes.reagent_size.height() / 2)
        second_reagent = AnswerField(second_reagent_point.x(), second_reagent_point.y())
        self._reagent_items.append(second_reagent)
        self.addToGroup(second_reagent)

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        for child_item in self.childItems():
            if child_item.sceneBoundingRect().contains(event.scenePos()):
                child_item.hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        for child_item in self.childItems():
            if child_item.sceneBoundingRect().contains(event.lastScenePos()):
                child_item.hoverLeaveEvent(event)

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        for child_item in self.childItems():
            if child_item.sceneBoundingRect().contains(event.scenePos()):
                child_item.mousePressEvent(event)

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        painter.save()

        self._paint_arrows(painter)

        painter.restore()

    def _paint_arrows(self, painter: QPainter) -> None:
        top_left_chain_point = self.boundingRect().topLeft()
        first_arrow_start_point = top_left_chain_point + QPointF(
            Sizes.add_item_size.width() + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap,
            Sizes.reagent_size.height() / 2
        )
        painter.drawEllipse(top_left_chain_point, 5, 5)

        for i in range(0, len(self._reagent_items) - 1):
            arrow_start_point = first_arrow_start_point + i * QPointF(
                Sizes.arrow_size.width() + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap,
                0
            )

            arrow_path = QPainterPath(arrow_start_point)
            arrow_next_point = arrow_start_point + QPointF(Sizes.arrow_size.width(), 0)
            arrow_path.lineTo(arrow_next_point)
            arrow_next_point += QPointF(-Sizes.arrow_size.width() / 2, -Sizes.arrow_size.height() / 2)
            arrow_path.lineTo(arrow_next_point)
            arrow_next_point += QPointF(Sizes.arrow_size.width() / 2, Sizes.arrow_size.height() / 2)
            arrow_path.lineTo(arrow_next_point)
            arrow_next_point += QPointF(-Sizes.arrow_size.width() / 2, Sizes.arrow_size.height() / 2)
            arrow_path.lineTo(arrow_next_point)

            painter.drawPath(arrow_path)
