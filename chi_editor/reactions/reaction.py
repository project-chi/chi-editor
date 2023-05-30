from PyQt6.QtWidgets import QGraphicsItemGroup, QGraphicsItem


class Reaction(QGraphicsItemGroup):

    def __init__(self, x: float, y: float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setPos(x, y)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFiltersChildEvents(False)
        self.setAcceptHoverEvents(False)

    def to_string(self):
        pass
