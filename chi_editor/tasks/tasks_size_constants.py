from typing import ClassVar

from PyQt6.QtCore import QSizeF, QPointF


class Sizes:
    add_item_size: ClassVar[QSizeF] = QSizeF(40, 40)
    reagent_size: ClassVar[QSizeF] = QSizeF(60, 60)
    reagent_boarder_width: ClassVar[float] = 1
    plus_size: ClassVar[QSizeF] = QSizeF(20, 20)
    arrow_size: ClassVar[QSizeF] = QSizeF(60, 30)
    default_gap: ClassVar[float] = 20
    side_items_offset: ClassVar[float] = arrow_size.width() / 2 + default_gap
    plus_right_offset: ClassVar[float] = side_items_offset + reagent_size.width() + default_gap
    plus_top_offset: ClassVar[float] = -1 * plus_size.height() / 2
    plus_top_point_right: ClassVar[QPointF] = QPointF(plus_right_offset, plus_top_offset)
    plus_left_offset: ClassVar[float] = -1 * (
            side_items_offset + reagent_size.width() + default_gap + plus_size.width())
    plus_top_point_left: ClassVar[QPointF] = QPointF(plus_left_offset, plus_top_offset)
