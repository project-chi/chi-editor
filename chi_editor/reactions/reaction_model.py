class ReactionModel:
    _reagents: list[str] = {}
    _products: list[str] = {}

    def add_reagent(self, reagent: str) -> None:
        self._reagents.append(reagent)

    def add_product(self, product: str) -> None:
        self._reagents.append(product)

    def pop_reagent(self) -> None:
        if len(self._reagents) == 0:
            pass

        self._reagents.pop()

    def pop_product(self) -> None:
        if len(self._products) == 0:
            pass

        self._products.pop()
