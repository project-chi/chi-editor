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
    QGraphicsRectItem,
    QGraphicsSceneMouseEvent
)

from chi_editor.reactions.reaction_model import ReactionModel
from chi_editor.reactions.reagent_adder import ReagentAdder, GrowthDirection
from chi_editor.reactions.size_constants import Sizes


class ReactionItem(QGraphicsItemGroup):
    _model: ReactionModel
    _reagent_items: list[QGraphicsItem] = []
    _product_items: list[QGraphicsItem] = []
    _add_reagent_item: QGraphicsEllipseItem
    _add_product_item: QGraphicsEllipseItem

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self._add_reagent_item = ReagentAdder(self._reagent_items, GrowthDirection.LEFT)
        self._add_reagent_item.setRect(QRectF(
            QPointF(-1 * (
                        Sizes.arrow_width / 2 + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap + Sizes.add_item_size.width()),
                    -1 * Sizes.add_item_size.height() / 2), Sizes.add_item_size))
        self._add_product_item = ReagentAdder(self._product_items, GrowthDirection.RIGHT)
        self._add_product_item.setRect(
            QRectF(QPointF(Sizes.arrow_width / 2 + Sizes.default_gap + Sizes.reagent_size.width() + Sizes.default_gap,
                           -1 * Sizes.add_item_size.height() / 2),
                   Sizes.add_item_size))
        self.setFiltersChildEvents(False)
        self.setAcceptHoverEvents(False)

        self.addToGroup(self._add_product_item)
        self.addToGroup(self._add_reagent_item)

        self._addInitialItems()

    def _addInitialItems(self) -> None:
        self._reagent_items.append(QGraphicsRectItem(QRectF(
            QPointF(-1 * (Sizes.side_items_offset + Sizes.reagent_size.width()), -1 * Sizes.reagent_size.height() / 2),
            Sizes.reagent_size)))
        self._addLastReagent()

        self._product_items.append(QGraphicsRectItem(
            QRectF(QPointF(Sizes.side_items_offset, -1 * Sizes.reagent_size.height() / 2), Sizes.reagent_size)))
        self._addLastProduct()

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        if self._add_reagent_item.sceneBoundingRect().contains(event.scenePos()):
            self._add_reagent_item.hoverEnterEvent(event)
        elif self._add_product_item.sceneBoundingRect().contains(event.scenePos()):
            self._add_product_item.hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        if self._add_reagent_item.sceneBoundingRect().contains(event.lastScenePos()):
            self._add_reagent_item.hoverLeaveEvent(event)
        elif self._add_product_item.sceneBoundingRect().contains(event.lastScenePos()):
            self._add_product_item.hoverLeaveEvent(event)

    def mousePressEvent(self, event: 'QGraphicsSceneMouseEvent') -> None:
        if self._add_reagent_item.sceneBoundingRect().contains(event.lastScenePos()):
            self._add_reagent_item.mousePressEvent(event)
            self._addLastReagent()
        elif self._add_product_item.sceneBoundingRect().contains(event.lastScenePos()):
            self._add_product_item.mousePressEvent(event)
            self._addLastProduct()

    def _addLastReagent(self) -> None:
        last_reagent = self._reagent_items.pop()
        self.addToGroup(last_reagent)
        self._reagent_items.append(last_reagent)

    def _addLastProduct(self) -> None:
        last_product = self._product_items.pop()
        self.addToGroup(last_product)
        self._product_items.append(last_product)

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        painter.save()

        self._paint_arrow(painter)

        painter.restore()

    def _paint_arrow(self, painter: QPainter) -> None:
        arrow_path = QPainterPath(QPointF(-1 * Sizes.arrow_width / 2, 0))
        arrow_path.lineTo(Sizes.arrow_width / 2, 0)
        arrow_path.lineTo(0, Sizes.arrow_height / 2)
        arrow_path.lineTo(Sizes.arrow_width / 2, 0)
        arrow_path.lineTo(0, -1 * Sizes.arrow_height / 2)

        painter.drawPath(arrow_path)
