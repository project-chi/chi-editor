from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsTextItem, QGraphicsItem

from ...bases.tool import Tool


class Text(Tool):
    def action(self, event: QGraphicsSceneMouseEvent) -> None:
        new_text = QGraphicsTextItem()
        new_text.setPos(event.scenePos())
        new_text.setPlainText('H')
        new_text.setFont(QFont('Impact'))
        new_text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)

        self.canvas.addItem(new_text)

    @property
    def asset(self) -> str:
        return 'text'

