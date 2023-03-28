from typing import TYPE_CHECKING

from PyQt6.QtGui import QAction, QIcon

from chi_editor.constants import RESOURCES

if TYPE_CHECKING:
    from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

    from chi_editor.canvas import Canvas


class Tool(QAction):
    def __init__(self, canvas: "Canvas", *args, **kwargs) -> "None":
        super().__init__(*args, **kwargs)
        self.canvas = canvas

        icon_path = RESOURCES / "assets" / "toolbar" / f"{self.picture}.png"
        if not icon_path.exists():
            raise ValueError(
                f"can't find assets for {type(self).__name__}, checked at "
                f"{icon_path}\n module reference: {self.__module__}"
            )

        self.setIcon(QIcon(str(icon_path)))
        self.setCheckable(True)

    def mouse_press_event(self, event: "QGraphicsSceneMouseEvent") -> "None":
        pass

    def mouse_move_event(self, event: "QGraphicsSceneMouseEvent") -> "None":
        pass

    def mouse_release_event(self, event: "QGraphicsSceneMouseEvent") -> "None":
        pass

    @property
    def picture(self) -> "str":
        """Returns the name of the icon for this tool without extension."""
        raise NotImplementedError
