from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsTextItem, QGraphicsItem, QGraphicsSceneHoverEvent

from ...bases.tool import Tool
from ...canvas import Canvas


class Text(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_text = TextItem()
        new_text.setPos(event.scenePos())

        self.canvas.addItem(new_text)

    @property
    def asset(self) -> str:
        return 'text'


class TextItem(QGraphicsTextItem):
    _canvas: Canvas

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setPlainText('H')
        self.setFont(QFont('Impact'))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setAcceptHoverEvents(True)

    def hoverEnterEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        HoverText(is_right=True, parent=self)
        HoverText(is_right=False, parent=self)

    def hoverLeaveEvent(self, event: 'QGraphicsSceneHoverEvent') -> None:
        for i in self.childItems():
            self.scene().removeItem(i)


class HoverText(QGraphicsTextItem):
    def __init__(self, *args, is_right: bool, parent: TextItem, **kwargs):
        super().__init__(*args, **kwargs)
        self.setPlainText('+')
        self.setFont(QFont('Impact'))
        self.setAcceptHoverEvents(True)
        self.setParentItem(parent)
        self.setPos(10 if is_right else -10, 0)
