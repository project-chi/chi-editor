from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom


class CreateCopper(CreateAtom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "Cu"

    @property
    def picture(self) -> str:
        return "copper"