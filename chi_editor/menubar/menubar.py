from PyQt6.QtWidgets import QMenuBar

from ..editor import Editor
from .menu_tools import menu_tools


class CanvasMenuBar(QMenuBar):
    _editor: Editor

    def __init__(self, *args, editor: Editor, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._editor = editor

        mode_menu = self.addMenu("Mode")
        for MenuTool in menu_tools:
            menu_tool = MenuTool(editor=editor, parent=mode_menu)
            mode_menu.addAction(menu_tool)
