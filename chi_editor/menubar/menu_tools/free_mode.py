from ...bases.menu_tool import MenuTool


class FreeMode(MenuTool):
    def exec(self) -> None:
        pass

    @property
    def asset(self) -> str:
        raise NotImplemented

    @property
    def text_asset(self) -> str:
        return "Enter free mode"
