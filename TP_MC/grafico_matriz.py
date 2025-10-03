import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

M  = np.loadtxt("matriz1.dat", dtype=int)

# Colores: -1 -> azul, 0 -> blanco, 1 -> rojo
cmap  = ListedColormap(["#1f77b4", "#ffffff", "#d62728"])
bounds = [-1.5, -0.5, 0.5, 1.5]
norm   = BoundaryNorm(bounds, cmap.N)

plt.figure(figsize=(4.5, 4.5))
im = plt.imshow(M, cmap=cmap, norm=norm, interpolation="nearest", aspect="equal")
cbar = plt.colorbar(im, ticks=[-1, 0, 1])
cbar.set_label("Valor")

plt.title("matriz.dat (2D)")
plt.xlabel("Columna")
plt.ylabel("Fila")

# (Opcional) cuadrícula suave entre celdas
for y in range(M.shape[0] + 1):
    plt.axhline(y - 0.5, color="k", lw=0.5, alpha=0.15)
for x in range(M.shape[1] + 1):
    plt.axvline(x - 0.5, color="k", lw=0.5, alpha=0.15)

plt.tight_layout()
plt.show()
