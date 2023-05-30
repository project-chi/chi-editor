from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom


class CreateZinc(CreateAtom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "Zn"

    @property
    def picture(self) -> str:
        return "zinc"
