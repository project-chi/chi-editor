from ...bases.menu_tool import MenuTool


class CreateTask(MenuTool):
    def exec(self) -> None:
        pass

    @property
    def asset(self) -> str:
        raise NotImplemented

    @property
    def text_asset(self) -> str:
        return "Create task"
