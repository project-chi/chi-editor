from typing import ClassVar

from PyQt6.QtCore import QSizeF, QPointF


class Sizes:
    add_item_size: ClassVar[QSizeF] = QSizeF(20, 20)
    reagent_size: ClassVar[QSizeF] = QSizeF(120, 120)
    plus_size: ClassVar[QSizeF] = QSizeF(20, 20)
    default_gap: ClassVar[float] = 20
    arrow_height: ClassVar[float] = 30
    arrow_width: ClassVar[float] = 60
    side_items_offset: ClassVar[float] = arrow_width / 2 + default_gap
    plus_right_offset: ClassVar[float] = side_items_offset + reagent_size.width() + default_gap
    plus_top_offset: ClassVar[float] = -1 * plus_size.height() / 2
    plus_top_point_right: ClassVar[QPointF] = QPointF(plus_right_offset, plus_top_offset)
    plus_left_offset: ClassVar[float] = -1 * (
                side_items_offset + reagent_size.width() + default_gap + plus_size.width())
    plus_top_point_left: ClassVar[QPointF] = QPointF(plus_left_offset, plus_top_offset)
