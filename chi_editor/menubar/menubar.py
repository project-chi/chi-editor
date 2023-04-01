from typing import TYPE_CHECKING

from PyQt6.QtWidgets import QMenuBar

from chi_editor.menubar.menu_tools import mode_tools, task_tools

if TYPE_CHECKING:
    from chi_editor.editor import Editor


class CanvasMenuBar(QMenuBar):
    _editor: "Editor"

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._editor = editor

        mode_menu = self.addMenu("Mode")
        for MenuTool in mode_tools:
            menu_tool = MenuTool(editor=editor, parent=mode_menu)
            mode_menu.addAction(menu_tool)

        task_menu = self.addMenu("Tasks")
        for MenuTool in task_tools:
            menu_tool = MenuTool(editor=editor, parent=task_menu)
            task_menu.addAction(menu_tool)
