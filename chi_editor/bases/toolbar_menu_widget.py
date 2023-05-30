from PyQt6.QtGui import QIcon, QAction
from PyQt6.QtWidgets import QToolButton

from chi_editor.bases.tool import Tool
from chi_editor.constants import RESOURCES


class ToolbarMenuWidget(QToolButton):
    _current_tool: QAction

    def __init__(self, *args, **kwargs) -> "None":
        super().__init__(*args, **kwargs)

        icon_path = RESOURCES / "assets" / "toolbar" / f"{self.picture}.png"
        if not icon_path.exists():
            raise ValueError(
                f"can't find assets for {type(self).__name__}, checked at "
                f"{icon_path}\n module reference: {self.__module__}"
            )

        self.setIcon(QIcon(str(icon_path)))
        self.triggered.connect(self._set_tool)
        self.setCheckable(True)

    def _set_tool(self, tool: Tool):
        self._current_tool = tool
        self.setIcon(tool.icon())

    def get_tool(self) -> QAction:
        return self._current_tool

    @property
    def picture(self) -> "str":
        """Returns the name of the icon for this menu without extension."""
        raise NotImplementedError
