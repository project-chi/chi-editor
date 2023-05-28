from typing import Optional

from PyQt6.QtCore import QRectF, QPointF
from PyQt6.QtGui import QPainter, QPainterPath, QColor
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QGraphicsItemGroup,
    QGraphicsEllipseItem,
    QGraphicsSceneHoverEvent,
    QStyleOptionGraphicsItem,
    QWidget,
    QGraphicsSceneMouseEvent
)

from chi_editor.reactions.reaction_model import ReactionModel
from chi_editor.reactions.reagent_adder import ReagentAdder, GrowthDirection
from chi_editor.reactions.size_constants import Sizes
from chi_editor.tasks.answer_field.answer_field import AnswerField


class ReactionItem(QGraphicsItemGroup):
    _model: ReactionModel
    _reagent_items: list[QGraphicsItem] = []
    _product_items: list[QGraphicsItem] = []
    _add_reagent_item: QGraphicsEllipseItem
    _add_product_item: QGraphicsEllipseItem

    def __init__(self, x: float, y: float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setPos(x, y)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFiltersChildEvents(False)
        self.setAcceptHoverEvents(False)

        self._addAddButtons()
        self._addInitialItems()

    def _addAddButtons(self) -> None:
        reagent_top_left_point = self.pos() + QPointF(-1 * (
                Sizes.arrow_width / 2 + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap + Sizes.add_item_size.width()),
                                                      -1 * Sizes.add_item_size.height() / 2)
        self._add_reagent_item = ReagentAdder(self._reagent_items, GrowthDirection.LEFT, reagent_top_left_point)

        product_top_left_point = self.pos() + QPointF(
            Sizes.arrow_width / 2 + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap,
            -1 * Sizes.add_item_size.height() / 2)
        self._add_product_item = ReagentAdder(self._product_items, GrowthDirection.RIGHT, product_top_left_point)

        self.addToGroup(self._add_product_item)
        self.addToGroup(self._add_reagent_item)

    def _addInitialItems(self) -> None:
        first_reagent_point = self.pos() + QPointF(-1 * (Sizes.side_items_offset + Sizes.reagent_size.width()),
                                                   -1 * Sizes.reagent_size.height() / 2)
        self._reagent_items.append(AnswerField(first_reagent_point.x(), first_reagent_point.y()))

        self._addLastReagent()

        first_product_point = self.pos() + QPointF(Sizes.side_items_offset, -1 * Sizes.reagent_size.height() / 2)
        self._product_items.append(AnswerField(first_product_point.x(), first_product_point.y()))
        self._addLastProduct()

    def _addLastReagent(self) -> None:
        last_reagent = self._reagent_items.pop()
        self.addToGroup(last_reagent)
        self._reagent_items.append(last_reagent)

    def _addLastProduct(self) -> None:
        last_product = self._product_items.pop()
        self.addToGroup(last_product)
        self._product_items.append(last_product)

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        print(self.sceneBoundingRect().contains(event.scenePos()))
        if self._add_reagent_item.sceneBoundingRect().contains(event.scenePos()):
            self._add_reagent_item.hoverEnterEvent(event)
        elif self._add_product_item.sceneBoundingRect().contains(event.scenePos()):
            self._add_product_item.hoverEnterEvent(event)

        for reagent_item in self._reagent_items:
            if reagent_item.sceneBoundingRect().contains(event.scenePos()):
                reagent_item.hoverEnterEvent(event)
        for product_item in self._product_items:
            if product_item.sceneBoundingRect().contains(event.scenePos()):
                product_item.hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        if self._add_reagent_item.sceneBoundingRect().contains(event.lastScenePos()):
            self._add_reagent_item.hoverLeaveEvent(event)
        elif self._add_product_item.sceneBoundingRect().contains(event.lastScenePos()):
            self._add_product_item.hoverLeaveEvent(event)

        for reagent_item in self._reagent_items:
            if reagent_item.sceneBoundingRect().contains(event.lastScenePos()):
                reagent_item.hoverLeaveEvent(event)
        for product_item in self._product_items:
            if product_item.sceneBoundingRect().contains(event.lastScenePos()):
                product_item.hoverLeaveEvent(event)

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        if self._add_reagent_item.sceneBoundingRect().contains(event.scenePos()):
            self._add_reagent_item.mousePressEvent(event)
            self._addLastReagent()
        elif self._add_product_item.sceneBoundingRect().contains(event.scenePos()):
            self._add_product_item.mousePressEvent(event)
            self._addLastProduct()

        for reagent_item in self._reagent_items:
            if reagent_item.sceneBoundingRect().contains(event.scenePos()):
                reagent_item.mousePressEvent(event)
        for product_item in self._product_items:
            if product_item.sceneBoundingRect().contains(event.scenePos()):
                product_item.mousePressEvent(event)

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        painter.save()

        self._paint_arrow(painter)
        self._paint_pluses(painter)

        painter.restore()

    def _paint_arrow(self, painter: QPainter) -> None:
        arrow_path = QPainterPath(QPointF(-1 * Sizes.arrow_width / 2, 0))
        arrow_path.lineTo(Sizes.arrow_width / 2, 0)
        arrow_path.lineTo(0, Sizes.arrow_height / 2)
        arrow_path.lineTo(Sizes.arrow_width / 2, 0)
        arrow_path.lineTo(0, -1 * Sizes.arrow_height / 2)

        painter.drawPath(arrow_path)

    def _paint_pluses(self, painter: QPainter) -> None:
        for i in range(0, len(self._reagent_items) - 1):
            plus_top_left = Sizes.plus_top_point_left - i * QPointF(
                Sizes.plus_size.width() + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap, 0)

            painter.drawLine(plus_top_left + QPointF(Sizes.plus_size.width() / 2, 0),
                             plus_top_left + QPointF(Sizes.plus_size.width() / 2, Sizes.plus_size.height()))

            painter.drawLine(plus_top_left + QPointF(0, Sizes.plus_size.height() / 2),
                             plus_top_left + QPointF(Sizes.plus_size.width(), Sizes.plus_size.height() / 2))

        for i in range(0, len(self._product_items) - 1):
            plus_top_left = Sizes.plus_top_point_right + i * QPointF(
                Sizes.plus_size.width() + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap, 0)

            painter.drawLine(plus_top_left + QPointF(Sizes.plus_size.width() / 2, 0),
                             plus_top_left + QPointF(Sizes.plus_size.width() / 2, Sizes.plus_size.height()))

            painter.drawLine(plus_top_left + QPointF(0, Sizes.plus_size.height() / 2),
                             plus_top_left + QPointF(Sizes.plus_size.width(), Sizes.plus_size.height() / 2))
