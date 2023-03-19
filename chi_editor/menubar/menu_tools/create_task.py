from ...bases.menu_tool import MenuTool
from ...editor_mode import EditorMode


class CreateTask(MenuTool):
    def exec(self) -> None:
        self._editor.setMode(EditorMode.CREATE_MODE)

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Create task"
