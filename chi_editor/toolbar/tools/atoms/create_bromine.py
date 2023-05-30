from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom


class CreateBromine(CreateAtom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "Br"

    @property
    def picture(self) -> str:
        return "bromine"
