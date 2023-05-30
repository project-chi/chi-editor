
from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom


class CreateCarbon(CreateAtom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "C"

    @property
    def picture(self) -> str:
        return "carbon"
