from ...bases.menu_tool import MenuTool


class SolveTask(MenuTool):
    def exec(self) -> None:
        pass

    @property
    def asset(self) -> str:
        raise NotImplemented

    @property
    def text_asset(self) -> str:
        return "Solve task"
