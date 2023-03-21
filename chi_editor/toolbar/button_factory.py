from enum import Enum
from pathlib import Path
from typing import Generator, Iterable

from PyQt6.QtGui import QAction, QIcon


class Tools(str, Enum):
    Arrow = "Arrow"
    Hand = "Hand"
    Bond = "Bond"
    Structure = "Structure"
    Atom = "Atom"
    Block = "Block"
    Reaction = "Reaction"
    Text = "Text"
    Eraser = "Eraser"


def create_buttons(query: Iterable[str]) -> Generator[QAction, None, None]:
    icon_library = Path(__file__).parent.parent.parent / "resources" / "toolbar"

    for action_name in query:
        icon_path = icon_library / f"{action_name.lower()}.png"
        if not icon_path.exists():
            raise ValueError(f"can\'t find assets for {action_name = }")

        action = QAction(QIcon(str(icon_path)), action_name)
        action.setCheckable(True)

        yield action
