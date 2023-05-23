from typing import Optional, ClassVar

from PyQt6.QtCore import QSizeF, QRectF, QPointF
from PyQt6.QtGui import QPainter, QPainterPath
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QGraphicsItemGroup,
    QGraphicsEllipseItem,
    QGraphicsSceneHoverEvent,
    QStyleOptionGraphicsItem,
    QWidget,
)

from chi_editor.reactions.reaction_model import ReactionModel
from chi_editor.reactions.reagent_adder import ReagentAdder, GrowthDirection


class ReactionItem(QGraphicsItemGroup):
    _model: ReactionModel
    _reagent_items: list[QGraphicsItem] = {}
    _product_items: list[QGraphicsItem] = {}
    _add_reagent_item: QGraphicsEllipseItem
    _add_product_item: QGraphicsEllipseItem
    add_item_size: ClassVar[QSizeF] = QSizeF(20, 20)
    reagent_size: ClassVar[QSizeF] = QSizeF(120, 120)
    plus_size: ClassVar[QSizeF] = QSizeF(20, 20)
    default_gap: ClassVar[float] = 20
    arrow_height: ClassVar[float] = 30
    arrow_width: ClassVar[float] = 60

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self._add_reagent_item = ReagentAdder(self._reagent_items, GrowthDirection.LEFT, parent=self)
        self._add_reagent_item.setRect(QRectF(
            QPointF(-1 * (self.arrow_width / 2 + self.default_gap + self.add_item_size.width()),
                    -1 * self.add_item_size.height() / 2), self.add_item_size))
        self._add_product_item = ReagentAdder(self._product_items, GrowthDirection.RIGHT, parent=self)
        self._add_product_item.setRect(
            QRectF(QPointF(self.arrow_width / 2 + self.default_gap, -1 * self.add_item_size.height() / 2),
                   self.add_item_size))
        self.setFiltersChildEvents(False)
        self.setAcceptHoverEvents(False)

        self.addToGroup(self._add_product_item)
        self.addToGroup(self._add_reagent_item)

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

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        painter.save()

        self._paint_arrow(painter)

        painter.restore()

    def _paint_arrow(self, painter: QPainter):
        arrow_path = QPainterPath(QPointF(-1 * self.arrow_width / 2, 0))
        arrow_path.lineTo(self.arrow_width / 2, 0)
        arrow_path.lineTo(0, self.arrow_height / 2)
        arrow_path.lineTo(self.arrow_width / 2, 0)
        arrow_path.lineTo(0, -1 * self.arrow_height / 2)

        painter.drawPath(arrow_path)
