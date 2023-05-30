from chi_editor.toolbar.tools.atoms.create_atom import CreateAtom
from chi_editor.toolbar.tools.atoms.create_hydrogen import CreateHydrogen
from chi_editor.toolbar.tools.atoms.create_helium import CreateHelium
from chi_editor.toolbar.tools.atoms.create_lithium import CreateLithium
from chi_editor.toolbar.tools.atoms.create_boron import CreateBoron
from chi_editor.toolbar.tools.atoms.create_carbon import CreateCarbon
from chi_editor.toolbar.tools.atoms.create_nitrogen import CreateNitrogen
from chi_editor.toolbar.tools.atoms.create_oxygen import CreateOxygen
from chi_editor.toolbar.tools.atoms.create_flourine import CreateFluorine
from chi_editor.toolbar.tools.atoms.create_sodium import CreateSodium
from chi_editor.toolbar.tools.atoms.create_magnesium import CreateMagnesium
from chi_editor.toolbar.tools.atoms.create_aluminium import CreateAluminium
from chi_editor.toolbar.tools.atoms.create_silicon import CreateSilicon
from chi_editor.toolbar.tools.atoms.create_phosphorous import CreatePhosphorous
from chi_editor.toolbar.tools.atoms.create_sulfur import CreateSulfur
from chi_editor.toolbar.tools.atoms.create_chlorine import CreateChlorine
from chi_editor.toolbar.tools.atoms.create_potassium import CreatePotassium
from chi_editor.toolbar.tools.atoms.create_calcium import CreateCalcium
from chi_editor.toolbar.tools.atoms.create_chromium import CreateChromium
from chi_editor.toolbar.tools.atoms.create_manganese import CreateManganese
from chi_editor.toolbar.tools.atoms.create_iron import CreateIron
from chi_editor.toolbar.tools.atoms.create_copper import CreateCopper
from chi_editor.toolbar.tools.atoms.create_zinc import CreateZinc
from chi_editor.toolbar.tools.atoms.create_bromine import CreateBromine

atom_tools: tuple[type[CreateAtom], ...] = (
    CreateHydrogen,
    CreateHelium,
    CreateLithium,
    CreateBoron,
    CreateCarbon,
    CreateNitrogen,
    CreateOxygen,
    CreateFluorine,
    CreateSodium,
    CreateMagnesium,
    CreateAluminium,
    CreateSilicon,
    CreatePhosphorous,
    CreateSulfur,
    CreateChlorine,
    CreatePotassium,
    CreateCalcium,
    CreateChromium,
    CreateManganese,
    CreateIron,
    CreateCopper,
    CreateZinc,
    CreateBromine
)
