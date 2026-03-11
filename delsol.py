# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 12:26:24 2026

@author: Александр
"""

from collections import defaultdict

input_pdb = "3_75_AG.pdb"
output_pdb = "3_75_AG_filtered2.pdb"

# диапазон по Z, в котором растворитель должен ОСТАТЬСЯ в ангстремах
z_min = 50.0
z_max = 200.0 # Для системы размером 4 z = 230

# имена воды и ионов
water_resnames = {"HOH", "WAT", "SOL", "TIP3"}
ion_resnames = {"NA", "CL", "K", "CA", "MG"}

ion_charge = {
    "NA": +1,
    "K":  +1,
    "CL": -1,
    "CA": +2,
    "MG": +2,
}

# ---------- читаем файл, делим на 3 части ----------
pre_atoms = []     # всё до первого ATOM
post_atoms = []    # всё после последнего ATOM
molecules = []     # [(resname, z, atoms)]

current = []
current_resname = None
current_resid = None
seen_atom = False

with open(input_pdb) as f:
    lines = f.readlines()

# найдём диапазон ATOM-блока
atom_indices = [i for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))]
first_atom = atom_indices[0]
last_atom = atom_indices[-1]

pre_atoms = lines[:first_atom]
post_atoms = lines[last_atom + 1:]

# ---------- парсим молекулы ----------
for line in lines[first_atom:last_atom + 1]:
    if not line.startswith(("ATOM", "HETATM")):
        continue

    resname = line[17:20].strip()
    resid = line[22:26].strip()

    if current_resname is None:
        current_resname = resname
        current_resid = resid

    if resname != current_resname or resid != current_resid:
        z = float(current[0][46:54])
        molecules.append((current_resname, z, current))
        current = []
        current_resname = resname
        current_resid = resid

    current.append(line)

if current:
    z = float(current[0][46:54])
    molecules.append((current_resname, z, current))

# ---------- удаление по Z ----------
kept = []
removed_water_mol = 0
removed_water_atoms = 0
removed_ions = defaultdict(int)
removed_charge = 0

for resname, z, atoms in molecules:

    if resname in water_resnames and (z < z_min or z > z_max):
        removed_water_mol += 1
        removed_water_atoms += len(atoms)
        continue

    if resname in ion_resnames and (z < z_min or z > z_max):
        removed_ions[resname] += 1
        removed_charge += ion_charge[resname]
        continue

    kept.append((resname, atoms))

# ---------- нейтрализация ----------
final_atoms = []
remaining_ions = []

for resname, atoms in kept:
    if resname in ion_resnames:
        remaining_ions.append((resname, atoms))
    else:
        final_atoms.append(atoms)

def sign(x):
    return 1 if x > 0 else -1

charge_to_fix = removed_charge

for resname, atoms in remaining_ions:
    q = ion_charge[resname]
    if charge_to_fix != 0 and sign(q) != sign(charge_to_fix):
        charge_to_fix += q
        removed_ions[resname] += 1
        continue
    final_atoms.append(atoms)

# ---------- запись (ВАЖНО) ----------
with open(output_pdb, "w") as out:
    for l in pre_atoms:
        out.write(l)
    for atoms in final_atoms:
        for l in atoms:
            out.write(l)
    for l in post_atoms:
        out.write(l)

# ---------- отчёт ----------
print("===================================")
print(f"Удалено воды : {removed_water_mol} молекул ({removed_water_atoms} атомов)")
print("Удалено ионов (всего):")
for ion, n in removed_ions.items():
    print(f"  {ion:>3} : {n}")
print("Финальный заряд системы: 0")
print("Результат:", output_pdb)
print("===================================")