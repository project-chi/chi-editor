from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom


class CreateManganese(CreateAtom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "Mn"

    @property
    def picture(self) -> str:
        return "manganese"
