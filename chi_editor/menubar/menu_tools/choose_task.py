from ...bases.menu_tool import MenuTool


class ChooseTask(MenuTool):
    def exec(self) -> None:
        self._editor.createChooseTaskDialog()

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Choose task"
