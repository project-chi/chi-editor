from PyQt6.QtGui import QAction, QIcon

from ..constants import RESOURCES
from ..editor import Editor


class MenuTool(QAction):
    _editor: Editor
    
    def __init__(self, *args, editor: Editor, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._editor = editor

        icon_path = RESOURCES / 'assets' / 'toolbar' / f'{self.asset}.png'
        if icon_path.exists():
            self.setIcon(QIcon(str(icon_path)))

        self.setText(self.text_asset)

        self.triggered().connect(self.exec())

    def exec(self) -> None:
        pass

    @property
    def asset(self) -> str:
        """Returns the name of the icon for this menu tool without extension."""
        raise NotImplemented

    @property
    def text_asset(self) -> str:
        """Returns the text of this tool."""
        raise NotImplemented
