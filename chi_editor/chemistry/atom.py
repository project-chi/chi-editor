class Atom:
    name: str
    x: float
    y: float

    def __init__(self, name: str, x: int, y: int) -> None:
        self.name = name
        self.x = x
        self.y = y
