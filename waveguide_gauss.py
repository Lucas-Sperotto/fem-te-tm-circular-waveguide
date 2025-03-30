import numpy as np
import matplotlib.pyplot as plt
import subprocess
import meshio
import pandas as pd
import os
from scipy.spatial import cKDTree
from scipy.linalg import eigh

# Gerar malha circular
def generate_circle_mesh(radius=0.01, filename='teste02'):
    geo_content = f"""
    lc = {radius/20};
    Point(1) = {{0, 0, 0, lc}};
    Point(2) = {{{radius}, 0, 0, lc}};
    Point(3) = {{0, {radius}, 0, lc}};
    Point(4) = {{-{radius}, 0, 0, lc}};
    Point(5) = {{0, -{radius}, 0, lc}};
    Circle(1) = {{2, 1, 3}};
    Circle(2) = {{3, 1, 4}};
    Circle(3) = {{4, 1, 5}};
    Circle(4) = {{5, 1, 2}};
    Line Loop(5) = {{1,2,3,4}};
    Plane Surface(6) = {{5}};
    Physical Curve(100) = {{1,2,3,4}};
    Physical Surface(200) = {{6}};
    Mesh 2;
    Save "{filename}.msh";
    """
    with open(f"{filename}.geo", 'w') as f:
        f.write(geo_content)
    subprocess.run(["gmsh", f"{filename}.geo", "-2", "-format", "msh2"], check=True)

# Funções peso para EFGM
def peso(r, dmax):
    q = r/dmax
    return np.where(q <= 1, (1 - q)**2, 0)

# Montagem das matrizes globais EFGM
def montar_matrizes_efgm(points, dmax):
    Nn = len(points)
    S = np.zeros((Nn, Nn))
    T = np.zeros((Nn, Nn))

    tree = cKDTree(points)

    for i, pt in enumerate(points):
        idx = tree.query_ball_point(pt, dmax)
        xi = points[idx]
        ri = np.linalg.norm(xi - pt, axis=1)
        w = peso(ri, dmax)

        # Matriz de momento
        P = np.column_stack((np.ones(len(idx)), xi[:, 0], xi[:, 1]))
        A = P.T @ (np.diag(w)) @ P
        A_inv = np.linalg.inv(A)

        # Derivadas das funções de base
        B = np.zeros((2, len(idx)))
        for j in range(len(idx)):
            dpdx = np.array([0, 1, 0])
            dpdy = np.array([0, 0, 1])
            phi_grad = A_inv @ (P[j] * w[j])
            B[:, j] = phi_grad[1:]

        # Montagem local
        S_local = B.T @ B
        T_local = np.eye(len(idx))  # Simplificação, considerar integração numérica

        # Montagem global
        for a, global_a in enumerate(idx):
            for b, global_b in enumerate(idx):
                S[global_a, global_b] += S_local[a, b]
                T[global_a, global_b] += T_local[a, b]

    return S, T

# Resolver EFGM
def solve_efgm(radius=0.01, dmax_factor=2.0, num_modes=24, filename='teste02'):
    generate_circle_mesh(radius, filename)
    mesh = meshio.read(f"{filename}.msh")
    points = mesh.points[:, :2]
    dmax = radius / dmax_factor

    S, T = montar_matrizes_efgm(points, dmax)

    vals, vecs = eigh(S, T)
    kc = np.sqrt(np.real(vals[:num_modes]))
    fc = 3e8 * kc / (2 * np.pi)

    output_dir = f'out/results/efgm_{filename}'
    os.makedirs(output_dir, exist_ok=True)
    df = pd.DataFrame({'Modo': np.arange(1, num_modes+1), 'Frequência (GHz)': fc/1e9})
    df.to_csv(os.path.join(output_dir, 'frequencias_modos.csv'), index=False)

    return fc

# Execução principal
if __name__ == "__main__":
    fc = solve_efgm()
    print("Frequências EFGM (GHz):", fc/1e9)