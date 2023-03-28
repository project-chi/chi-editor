from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFocusEvent, QFont, QKeyEvent
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QGraphicsSceneHoverEvent,
    QGraphicsSceneMouseEvent,
    QGraphicsTextItem,
)

from ...bases.tool import Tool


class Text(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_text = TextItem()
        new_text.setPos(event.scenePos())

        self.canvas.addItem(new_text)

    @property
    def picture(self) -> str:
        return "text"


class TextItem(QGraphicsTextItem):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setPlainText("H")
        self.setFont(QFont("Impact"))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setAcceptHoverEvents(True)

    def hoverEnterEvent(self, event: QGraphicsSceneHoverEvent) -> None:
        HoverText(is_right=True, parent=self)
        HoverText(is_right=False, parent=self)

    def hoverLeaveEvent(self, event: QGraphicsSceneHoverEvent) -> None:
        for i in self.childItems():
            self.scene().removeItem(i)

    def mouseDoubleClickEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditable)
        self.setFocus(Qt.FocusReason.MouseFocusReason)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if event.key() != Qt.Key.Key_Return:
            QGraphicsTextItem.keyPressEvent(self, event)
        else:
            self.clearFocus()
            self.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)

    def focusOutEvent(self, event: QFocusEvent) -> None:
        self.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)


class HoverText(QGraphicsTextItem):
    def __init__(self, *args, is_right: bool, parent: TextItem, **kwargs):
        super().__init__(*args, **kwargs)
        self.setPlainText("+")
        self.setFont(QFont("Impact"))
        self.setAcceptHoverEvents(True)
        self.setParentItem(parent)
        self.setPos(parent.boundingRect().width() - 7 if is_right else -10, 0)
