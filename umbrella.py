# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 14:22:40 2026

@author: Александр
"""
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# =========================
# ПАРАМЕТРЫ
# =========================

pdb_file = "INH_3_75_AG.pdb"

energy_files = [
    "G_0_10.xvg",
    "G_5_10.xvg",
    "G_10_15.xvg",
    "G_10_20.xvg",
    "G_15_20.xvg",
    "G_20_25.xvg",
    "G_20_30.xvg",
    "G_25_30.xvg"
]

energy_labels = [
    "G_0_10",
    "G_5_10",
    "G_10_15",
    "G_10_20",
    "G_15_20",
    "G_20_25",
    "G_20_30",
    "G_25_30"
]

energy_colors = [
    "black",
    "tab:orange",
    "tab:green",
    "tab:blue",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:cyan"
]


z_center_abs = 7.39   # <-- задаёшь сам (nm)

nbins = 300
threshold_fraction = 0.05

components = {
    "MurNGc-GlcNAc": "resname GNA or resname MNG or resname MNGA or resname MNGB or resname MNGP",
    "Peptide Bridges": "resname PBDD or resname PBDA",
    "Arabinogalactan": "resname PGAG",
    "Water": "resname SOL"
}

colors = {
    "MurNGc-GlcNAc": "blue",
    "Peptide Bridges": "green",
    "Arabinogalactan": "red",
    "Water": "gray"
}


amu_to_kg = 1.66053906660e-27
nm3_to_m3 = 1e-27


# =========================
# РАСЧЁТ ПЛОТНОСТЕЙ (ДЛЯ ОБЛАСТЕЙ)
# =========================

u = mda.Universe(pdb_file)

Lx_nm = u.dimensions[0] * 0.1
Ly_nm = u.dimensions[1] * 0.1
area_nm2 = Lx_nm * Ly_nm

z_all_nm = u.atoms.positions[:, 2] * 0.1
z_min = z_all_nm.min()
z_max = z_all_nm.max()

bins = np.linspace(z_min, z_max, nbins + 1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
dz = bins[1] - bins[0]

densities = {}

for name, sel in components.items():
    atoms = u.select_atoms(sel)
    if len(atoms) == 0:
        print(f"{name}: отсутствует")
        continue

    z_nm = atoms.positions[:, 2] * 0.1
    masses = atoms.masses

    hist_mass_amu, _ = np.histogram(z_nm, bins=bins, weights=masses)
    hist_mass_kg = hist_mass_amu * amu_to_kg
    bin_volume_m3 = area_nm2 * dz * nm3_to_m3

    density = hist_mass_kg / bin_volume_m3
    densities[name] = density

# =========================
# ПЕРЕВОД В ОТНОСИТЕЛЬНЫЕ Z
# =========================

z_rel = bin_centers - z_center_abs


# =========================
# ОПРЕДЕЛЕНИЕ ОБЛАСТЕЙ
# =========================

regions = {}

for name, density in densities.items():
    threshold = threshold_fraction * np.max(density)
    mask = density > threshold
    regions[name] = mask


# =========================
# ЧТЕНИЕ XVG
# =========================

def read_xvg(filename):
    data = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            data.append([float(x) for x in line.split()])
    return np.array(data)


# =========================
# ПОСТРОЕНИЕ
# =========================

fig, ax = plt.subplots(figsize=(8, 5))

# =========================
# ОБЛАСТИ
# =========================

region_handles = []

for name, density in densities.items():

    threshold = threshold_fraction * np.max(density)
    mask = density > threshold

    if not np.any(mask):
        continue

    indices = np.where(mask)[0]
    segments = np.split(indices, np.where(np.diff(indices) != 1)[0] + 1)

    for segment in segments:
        z_start = z_rel[segment[0]]
        z_end   = z_rel[segment[-1]]

        ax.axvspan(
            z_start,
            z_end,
            color=colors.get(name, None),
            alpha=0.10
        )

    # создаём элемент легенды для области
    region_handles.append(
        Patch(
            facecolor=colors.get(name, None),
            alpha=0.3,
            label=name
        )
    )

# =========================
# ЭНЕРГИИ
# =========================



energy_handles = []

for file, label, color in zip(energy_files, energy_labels, energy_colors):
    data = read_xvg(file)
    z_energy = data[:, 0]
    energy = data[:, 1]

    # --- извлекаем числа из имени файла ---
    name = file.replace(".xvg", "")
    _, a, b = name.split("_")

    a = int(a)
    b = int(b)

    # --- выбираем стиль линии ---
    linestyle = "--" if abs(a - b) == 10 else "-"

    line, = ax.plot(
        z_energy,
        energy,
        label=label,
        color=color,
        linestyle=linestyle,
        linewidth=2
    )

    energy_handles.append(line)

ax.set_xlabel("z relative to center (nm)")
ax.set_ylabel("Energy (kJ/mol)")
ax.grid(alpha=0.3)

# =========================
# ДВЕ ЛЕГЕНДЫ
# =========================

legend_regions = ax.legend(
    handles=region_handles,
    loc="lower right",
    frameon=True
)

legend_energy = ax.legend(
    handles=energy_handles,
    loc="lower left",
    frameon=True
)

# важно: вернуть первую легенду
ax.add_artist(legend_regions)

plt.tight_layout()
plt.savefig("INH_3_AG_umb.png", dpi=300)
plt.show()