from ...bases.menu_tool import MenuTool

from ...editor import Editor


class SolveTask(MenuTool):
    def exec(self) -> None:
        self._editor.setMode(Editor.EditorMode.SOLVE_MODE)

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Solve task"
