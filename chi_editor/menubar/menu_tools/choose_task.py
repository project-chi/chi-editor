from ...bases.menu_tool import MenuTool


class ChooseTask(MenuTool):
    def exec(self) -> None:
        self._editor.openChooseTaskDialog()

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Choose task"
