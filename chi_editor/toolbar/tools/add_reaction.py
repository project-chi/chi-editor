from PyQt6.QtWidgets import (
    QGraphicsSceneMouseEvent,
)

from ...bases.tool import Tool
from chi_editor.reactions.reaction_item_group import ReactionItemGroup


class AddReaction(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_reaction = ReactionItemGroup(event.scenePos().x(), event.scenePos().y())

        self.canvas.addItem(new_reaction)

    @property
    def picture(self) -> str:
        return "reaction"
