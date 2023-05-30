from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom


class CreatePhosphorous(CreateAtom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "P"

    @property
    def picture(self) -> str:
        return "phosphorous"
