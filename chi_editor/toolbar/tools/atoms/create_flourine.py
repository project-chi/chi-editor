from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom


class CreateFluorine(CreateAtom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "F"

    @property
    def picture(self) -> str:
        return "fluorine"
