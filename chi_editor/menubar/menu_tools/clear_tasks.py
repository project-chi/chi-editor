from chi_editor.bases.menu_tool import MenuTool


class ClearTasks(MenuTool):
    def exec(self) -> None:
        # self._editor.server.deleteAllTasks()
        pass

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Clear tasks"
