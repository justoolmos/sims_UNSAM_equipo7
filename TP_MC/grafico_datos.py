import pandas as pd
import matplotlib.pyplot as plt

fname = "out.dat"

# Intento 1: asumir que el archivo tiene encabezado correcto
try:
    df = pd.read_csv(fname, engine="python", sep=r"\s*,\s*", comment="#")
    if not {"Paso","E_tot"}.issubset(df.columns):
        raise ValueError("Encabezado inesperado")
except Exception:
    # Intento 2: sin encabezado, asignar nombres manualmente
    df = pd.read_csv(fname, engine="python", sep=r"\s*,\s*", header=None)
    df.columns = ["Paso","E_tot","E_med","E_sqr","M_tot","M_med","M_sqr"]

# Asegurar tipos numéricos
df["Paso"]  = pd.to_numeric(df["Paso"], errors="coerce")
df["E_tot"] = pd.to_numeric(df["E_tot"], errors="coerce")
df = df.dropna(subset=["Paso","E_tot"])

# Graficar
plt.figure(figsize=(6,4))
plt.plot(df["Paso"], df["E_tot"], marker="o")
plt.xlabel("Paso")
plt.ylabel("E_tot")
plt.title("E_tot vs Paso")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# (Opcional) Guardar:
# plt.savefig("E_tot_vs_Paso.png", dpi=300, bbox_inches="tight")

