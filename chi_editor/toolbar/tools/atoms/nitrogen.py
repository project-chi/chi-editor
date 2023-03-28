
from ..atom import Atom


class Nitrogen(Atom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "N"

    @property
    def picture(self) -> str:
        return "nitrogen"
