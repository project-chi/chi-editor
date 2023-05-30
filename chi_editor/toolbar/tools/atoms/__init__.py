from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom
from chi_editor.toolbar.tools.atoms.create_carbon import CreateCarbon
from chi_editor.toolbar.tools.atoms.create_nitrogen import CreateNitrogen
from chi_editor.toolbar.tools.atoms.create_oxygen import CreateOxygen

atom_tools: tuple[type[CreateAtom], ...] = (
    CreateCarbon,
    CreateNitrogen,
    CreateOxygen
)
