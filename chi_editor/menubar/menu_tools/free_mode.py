from ...bases.menu_tool import MenuTool

from ...editor import Editor


class FreeMode(MenuTool):
    def exec(self) -> None:
        self._editor.setMode(Editor.EditorMode.FREE_MODE)

    @property
    def asset(self) -> str:
        raise NotImplemented

    @property
    def text_asset(self) -> str:
        return "Enter free mode"
