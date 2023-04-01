from chi_editor.bases.menu_tool import MenuTool
from chi_editor.editor_mode import EditorMode


class SolveTask(MenuTool):
    def exec(self) -> None:
        self._editor.setMode(EditorMode.SOLVE_MODE)

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Solve task"
