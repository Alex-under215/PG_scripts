# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 10:43:17 2026

@author: Александр
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------- параметры ----------
pdb_file = "INH_3_75_AG.pdb"

# selections MDAnalysis
components = {
    "Peptidoglycan": "not water and not (resname NA or resname CL)",
    "MurNGc-GlcNAc": "resname GNA or resname MNG or resname MNGA or resname MNGB or resname MNGP",
    "Peptide Bridges": "resname PBDD or resname PBDA",
    "Arabinogalactan": "resname PGAG",
    "Water": "resname SOL"
}

colors = {
    "Peptidoglycan": "blue",
    "MurNGc-GlcNAc": "orange",
    "Peptide Bridges": "green",
    "Arabinogalactan": "red",
    "Water": "purple"
}

nbins = 200

# ---------- константы ----------
amu_to_kg = 1.66053906660e-27   # kg
nm3_to_m3 = 1e-27              # m³

# ---------- папка со скриптом ----------
script_dir = Path(__file__).resolve().parent

# ищем все pdb в папке
pdb_files = sorted(script_dir.glob("*.pdb"))

if not pdb_files:
    print("❌ В папке нет PDB файлов")
    exit()

print(f"Найдено PDB файлов: {len(pdb_files)}")

# ---------- цикл по всем pdb ----------
for pdb_file in pdb_files:
    print(f"Обрабатываю: {pdb_file.name}")

    u = mda.Universe(pdb_file)

    # размеры ячейки (Å → nm)
    if u.dimensions is None:
        print(f"⚠ В {pdb_file.name} нет размеров ячейки — пропуск")
        continue

    Lx_nm = u.dimensions[0] * 0.1
    Ly_nm = u.dimensions[1] * 0.1
    area_nm2 = Lx_nm * Ly_nm

    # координаты Z (Å → nm)
    z_all_nm = u.atoms.positions[:, 2] * 0.1
    z_min = z_all_nm.min()
    z_max = z_all_nm.max()

    bins = np.linspace(z_min, z_max, nbins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    dz = bins[1] - bins[0]   # nm

    plt.figure(figsize=(8, 5))

    for name, sel in components.items():
        atoms = u.select_atoms(sel)
        
        # если группы нет — не рисуем
        if len(atoms) == 0:
            print(f"   - группа '{name}' отсутствует, пропуск")
            continue

        z_nm = atoms.positions[:, 2] * 0.1
        masses_amu = atoms.masses

        hist_mass_amu, _ = np.histogram(z_nm, bins=bins, weights=masses_amu)
        hist_mass_kg = hist_mass_amu * amu_to_kg

        bin_volume_m3 = area_nm2 * dz * nm3_to_m3
        density = hist_mass_kg / bin_volume_m3

        plt.plot(bin_centers, density, label=name, color=colors.get(name, None))

    # ---------- оформление ----------
    plt.xlabel("Z (nm)")
    plt.ylabel("Mass density (kg/m³)")
    plt.title(pdb_file.stem)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # сохраняем график
    output_png = script_dir / f"{pdb_file.stem}_density.png"
    plt.savefig(output_png, dpi=300)
    plt.close()

    print(f"  → сохранён график: {output_png.name}")

print("✅ Готово.")

atoms = u.select_atoms("resname MNGA or resname MNGB")
print("Atoms:", len(atoms))
print("Residues:", len(atoms.residues))