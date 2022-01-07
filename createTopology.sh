cat PSS.atomtypes.top >> solvated.top
cat PSS.moleculetype >> solvated.top
cat PSS.atoms.top >> solvated.top
cat PSS.bonds.top >> solvated.top
cat PSS.angles.top >> solvated.top
cat PSS.dihedrals.top >> solvated.top
head BUTA.atomtypes.top >> solvated.top
head BUTA.moleculetype >> solvated.top
head -9 BUTA.atoms.top >> solvated.top
head -9 BUTA.bonds.top >> solvated.top
head -10 BUTA.angles.top >> solvated.top
head -5 BUTA.dihedrals.top >> solvated.top
head WAT.atomtypes.top >> solvated.top
head WAT.moleculetype >> solvated.top
head -4 WAT.atoms.top >> solvated.top
head -3 WAT.bonds.top >> solvated.top
head -2 WAT.angles.top >> solvated.top
head WAT.dihedrals.top >> solvated.top
