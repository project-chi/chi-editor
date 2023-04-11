from chi_editor.bases.menu_tool import MenuTool

from chi_editor.dialog_windows.choose_task.local_task_dialog import LocalTaskDialog


class ChooseLocalTask(MenuTool):
    _choose_dialog: LocalTaskDialog

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._choose_dialog = LocalTaskDialog(editor=self._editor)
        self._choose_dialog.setModal(True)

    def exec(self) -> None:
        self._choose_dialog.exec()

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Choose task locally"
