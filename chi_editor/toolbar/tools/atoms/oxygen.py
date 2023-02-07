from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ..atom import Atom


class Oxygen(Atom):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._element = 'O'

    @property
    def asset(self) -> str:
        return 'oxygen'
