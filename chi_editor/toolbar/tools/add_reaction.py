from PyQt6.QtWidgets import (
    QGraphicsSceneMouseEvent,
)

from ...bases.tool import Tool
from chi_editor.reactions.reaction_item import ReactionItem


class AddReaction(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_reaction = ReactionItem()
        new_reaction.setPos(event.scenePos())

        self.canvas.addItem(new_reaction)

    @property
    def picture(self) -> str:
        return "text"
