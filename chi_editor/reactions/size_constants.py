from typing import ClassVar

from PyQt6.QtCore import QSizeF


class Sizes:
    add_item_size: ClassVar[QSizeF] = QSizeF(20, 20)
    reagent_size: ClassVar[QSizeF] = QSizeF(120, 120)
    plus_size: ClassVar[QSizeF] = QSizeF(20, 20)
    default_gap: ClassVar[float] = 20
    arrow_height: ClassVar[float] = 30
    arrow_width: ClassVar[float] = 60
    side_items_offset: ClassVar[float] = arrow_width / 2 + default_gap
