from chi_editor.toolbar.tools.bonds.create_bond import CreateBond
from chi_editor.toolbar.tools.bonds.create_single_bond import CreateSingleBond
from chi_editor.toolbar.tools.bonds.create_double_bond import CreateDoubleBond
from chi_editor.toolbar.tools.bonds.create_triple_bond import CreateTripleBond

bond_tools: tuple[type[CreateBond], ...] = (
    CreateSingleBond,
    CreateDoubleBond,
    CreateTripleBond,
)
