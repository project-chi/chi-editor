from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import QToolButton

from chi_editor.bases.tool import Tool


class ToolbarMenuWidget(QToolButton):
    _current_tool: QAction

    def __init__(self, *args, **kwargs) -> "None":
        super().__init__(*args, **kwargs)
        self.triggered.connect(self._set_tool)

    def _set_tool(self, tool: Tool):
        self._current_tool = tool
        self.setIcon(tool.icon())

    def get_tool(self) -> QAction:
        return self._current_tool