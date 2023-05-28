from PyQt6.QtWidgets import (
    QGraphicsSceneMouseEvent,
)

from ...bases.tool import Tool
from chi_editor.chains.chain_item_group import ChainItemGroup


class AddChain(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_chain = ChainItemGroup(event.scenePos().x(), event.scenePos().y())

        self.canvas.addItem(new_chain)

    @property
    def picture(self) -> str:
        return "text"
