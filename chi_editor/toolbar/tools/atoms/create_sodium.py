from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom


class CreateSodium(CreateAtom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "Na"

    @property
    def picture(self) -> str:
        return "sodium"
