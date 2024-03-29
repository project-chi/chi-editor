
from ..atom import Atom


class Carbon(Atom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = "C"

    @property
    def picture(self) -> str:
        return "carbon"
