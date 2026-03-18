# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 09:32:39 2026

@author: Александр
"""

#!/usr/bin/env python3
import numpy as np
from collections import defaultdict

# ===================== НАСТРОЙКИ =====================

SYSTEM_PDB = "4_75_AG_center1.pdb"      # куда вставляем
LIGAND_PDB = "INH.pdb"                  # что вставляем
OUTPUT_PDB = "INH_4_75_AG1.pdb"         # финальный файл

CUTOFF = 3.0  # Å

# координаты, куда вставлять (Å!)
Z_TARGET = 165.0  # Å — нужная координата Z

WATER_RESNAMES = {"SOL", "HOH", "TIP3", "WAT"}
REQUIRED_ATOMS = {"OW", "HW1", "HW2"}

ATOMIC_MASS = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "P": 30.974,
    "S": 32.06,
}

# ===================== БАЗОВЫЕ УТИЛИТЫ =====================

def parse_coord(line):
    return np.array([
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54])
    ])

def format_coord(line, coord):
    return (
        line[:30]
        + f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
        + line[54:]
    )

def get_resname(line):
    return line[17:20].strip()

def get_element(line):
    return line[76:78].strip()

# ===================== COM СИСТЕМЫ =====================

def center_of_mass_system(pdb):
    coords = []
    with open(pdb) as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")) and get_resname(line) not in WATER_RESNAMES:
                try:
                    coords.append(parse_coord(line))
                except ValueError:
                    continue
    if not coords:
        raise RuntimeError("Нет атомов для вычисления COM системы")
    return np.mean(np.array(coords), axis=0)

# ===================== ЛИГАНД =====================

def read_ligand(pdb):
    lines = []
    coords = []

    with open(pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                lines.append(line.rstrip())
                coords.append(parse_coord(line))

    return lines, np.array(coords)

def shift_ligand(lines, coords, target_xy, target_z):
    com = np.mean(coords, axis=0)

    shift = np.array([
        target_xy[0] - com[0],
        target_xy[1] - com[1],
        target_z     - com[2]
    ])

    new_lines = []
    new_coords = []

    for i, line in enumerate(lines):
        c = coords[i] + shift
        new_lines.append(format_coord(line, c))
        new_coords.append(c)

    return new_lines, np.array(new_coords)

# ===================== СБОР ВОДЫ (БЕЗ ID!) =====================

def collect_water_molecules(pdb):
    waters = []
    cur_coords = []
    cur_lines = []

    with open(pdb) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if get_resname(line) not in WATER_RESNAMES:
                continue

            atom = line[12:16].strip()

            if atom == "OW":
                if cur_coords:
                    waters.append({
                        "coords": np.array(cur_coords),
                        "lines": cur_lines
                    })
                cur_coords = []
                cur_lines = []

            cur_coords.append(parse_coord(line))
            cur_lines.append(line.rstrip())

    if cur_coords:
        waters.append({
            "coords": np.array(cur_coords),
            "lines": cur_lines
        })

    return waters

# ===================== УДАЛЕНИЕ ВОДЫ ВОКРУГ ЛИГАНДА =====================

def waters_to_remove(waters, ligand_coords, cutoff):
    remove_ids = set()

    for i, w in enumerate(waters):
        d = np.linalg.norm(
            w["coords"][:, None, :] - ligand_coords[None, :, :],
            axis=2
        )
        if np.any(d < cutoff):
            remove_ids.add(i)

    return remove_ids

# ===================== БЛИЖАЙШИЙ АТОМ К COM =====================

def nearest_atom_to_com(pdb, com):
    min_d = np.inf
    closest = None

    with open(pdb) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if get_resname(line) in WATER_RESNAMES:
                continue
            if get_element(line) == "H":
                continue

            d = np.linalg.norm(parse_coord(line) - com)
            if d < min_d:
                min_d = d
                closest = line.rstrip()

    return min_d, closest

# ===================== ЗАПИСЬ ФИНАЛЬНОГО PDB =====================

def write_pdb(system_pdb, output_pdb, waters, remove_ids, ligand_lines):
    remove_lines = set()
    for i in remove_ids:
        for l in waters[i]["lines"]:
            remove_lines.add(l)

    with open(system_pdb) as sys, open(output_pdb, "w") as out:
        for line in sys:
            if line.rstrip() in remove_lines:
                continue
            if line.startswith(("TER", "END", "ENDMDL")):
                continue
            out.write(line)

        for l in ligand_lines:
            out.write(l + "\n")

        out.write("END\n")

# ===================== MAIN =====================

def main():
    com = center_of_mass_system(SYSTEM_PDB)
    print("COM системы (без воды):", com)

    d, atom = nearest_atom_to_com(SYSTEM_PDB, com)
    print(f"Ближайший неводородный атом: {d:.3f} Å")
    print(atom)

    lig_lines, lig_coords = read_ligand(LIGAND_PDB)
    lig_lines, lig_coords = shift_ligand(
        lig_lines, lig_coords, com[:2], Z_TARGET
    )

    waters = collect_water_molecules(SYSTEM_PDB)
    remove_ids = waters_to_remove(waters, lig_coords, CUTOFF)

    print("Удалено молекул воды:", len(remove_ids))
    print("Удалено строк воды:",
          sum(len(waters[i]["lines"]) for i in remove_ids))

    write_pdb(SYSTEM_PDB, OUTPUT_PDB, waters, remove_ids, lig_lines)
    print("Файл записан:", OUTPUT_PDB)

if __name__ == "__main__":
    main()