from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QKeyEvent, QFocusEvent
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsTextItem, QGraphicsItem

from ...bases.tool import Tool

class MyText(QGraphicsTextItem):

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

class Text(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_text = MyText()
        new_text.setPos(event.scenePos())
        new_text.setPlainText('H')
        new_text.setFont(QFont('Impact'))
        new_text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)

        self.canvas.addItem(new_text)

    @property
    def asset(self) -> str:
        return 'text'