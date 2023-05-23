from PyQt6.QtCore import QSizeF, QRectF, QPointF
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QGraphicsItemGroup,
    QGraphicsEllipseItem,
    QGraphicsSceneHoverEvent
)

from chi_editor.reactions.reaction_model import ReactionModel
from chi_editor.reactions.reagent_adder import ReagentAdder

_default_add_item_size = QSizeF(40, 40)


class ReactionItem(QGraphicsItemGroup):
    _model: ReactionModel
    _reagent_items: list[QGraphicsItem] = {}
    _product_items: list[QGraphicsItem] = {}
    _add_reagent_item: QGraphicsEllipseItem
    _add_product_item: QGraphicsEllipseItem

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self._add_reagent_item = ReagentAdder(self._reagent_items, parent=self)
        self._add_reagent_item.setRect(QRectF(QPointF(-100, 0), _default_add_item_size))
        self._add_product_item = ReagentAdder(self._product_items, parent=self)
        self._add_product_item.setRect(QRectF(QPointF(100, 0), _default_add_item_size))
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
