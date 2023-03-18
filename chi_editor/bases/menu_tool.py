from PyQt6.QtGui import QAction, QIcon

from ..canvas import Canvas
from ..constants import RESOURCES


class MenuTool(QAction):
    def __init__(self, *args, canvas: Canvas, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.canvas = canvas

        icon_path = RESOURCES / 'assets' / 'toolbar' / f'{self.asset}.png'
        if icon_path.exists():
            self.setIcon(QIcon(str(icon_path)))

    def activate(self, event: QAction.ActionEvent) -> None:
        pass

    @property
    def asset(self) -> str:
        """Returns the name of the icon for this tool without extension."""
        raise NotImplemented
