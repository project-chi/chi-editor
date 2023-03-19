from ...bases.menu_tool import MenuTool

from ...editor import Editor


class CreateTask(MenuTool):
    def exec(self) -> None:
        self._editor.setMode(Editor.EditorMode.CREATE_MODE)

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Create task"
