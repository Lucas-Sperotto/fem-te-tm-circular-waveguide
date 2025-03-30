import csv
from scipy.special import jn_zeros, jnp_zeros

# Número de ordens m e raízes n
max_m = 10
num_zeros = 10

# Arquivos de saída
csv_tm_path = "out/bessel_roots_tm.csv"
csv_te_path = "out/bessel_roots_te.csv"

with open(csv_tm_path, mode="w", newline="") as file_tm, open(csv_te_path, mode="w", newline="") as file_te:
    writer_tm = csv.writer(file_tm)
    writer_te = csv.writer(file_te)

    # Cabeçalhos
    writer_tm.writerow(["m", "n", "Raiz"])
    writer_te.writerow(["m", "n", "Raiz"])

    for m in range(0, max_m + 1):
        # TM (Jm(x) = 0)
        zeros_TM = jn_zeros(m, num_zeros)
        for n, root in enumerate(zeros_TM, start=1):
            writer_tm.writerow([m, n, f"{root:.15f}"])

        # TE (Jm'(x) = 0)
        zeros_TE = jnp_zeros(m, num_zeros)
        for n, root in enumerate(zeros_TE, start=1):
            writer_te.writerow([m, n, f"{root:.15f}"])

print("Arquivos 'bessel_roots_tm.csv' e 'bessel_roots_te.csv' gerados com sucesso!")
