import csv
from scipy.special import jn_zeros, jnp_zeros

# Número de ordens m e raízes n
max_m = 10
num_zeros = 10

# Listas temporárias para armazenar os dados
tm_data = []
te_data = []

for m in range(0, max_m + 1):
    # TM (Jm(x) = 0)
    zeros_TM = jn_zeros(m, num_zeros)
    for n, root in enumerate(zeros_TM, start=1):
        tm_data.append((m, n, root))

    # TE (Jm'(x) = 0)
    zeros_TE = jnp_zeros(m, num_zeros)
    for n, root in enumerate(zeros_TE, start=1):
        te_data.append((m, n, root))

# Ordenar por raiz
tm_data.sort(key=lambda x: x[2])
te_data.sort(key=lambda x: x[2])

# Arquivos de saída
csv_tm_path = "bessel_roots_tm.csv"
csv_te_path = "bessel_roots_te.csv"

with open(csv_tm_path, mode="w", newline="") as file_tm:
    writer_tm = csv.writer(file_tm)
    writer_tm.writerow(["m", "n", "Raiz"])
    for row in tm_data:
        writer_tm.writerow([row[0], row[1], f"{row[2]:.15f}"])

with open(csv_te_path, mode="w", newline="") as file_te:
    writer_te = csv.writer(file_te)
    writer_te.writerow(["m", "n", "Raiz"])
    for row in te_data:
        writer_te.writerow([row[0], row[1], f"{row[2]:.15f}"])

print("Arquivos ordenados 'bessel_roots_tm.csv' e 'bessel_roots_te.csv' gerados com sucesso!")
