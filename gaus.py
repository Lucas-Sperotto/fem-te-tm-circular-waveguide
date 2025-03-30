import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
import meshio
import subprocess

# Função para gerar malha circular

def generate_circle_mesh(radius=0.01, filename='circle'):
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

# Função para montar matrizes locais com integração de Gauss (1 ponto)

def montar_matrizes_gauss_1pt(tri, points):
    # Coordenadas locais do ponto de Gauss no triângulo de referência
    xi, eta, w = 1/3, 1/3, 0.5  # peso = área do triângulo referência (0.5)

    # Derivadas das funções de base no espaço referência
    dphi_ref = np.array([[-1, -1], [1, 0], [0, 1]])

    # Coordenadas globais dos vértices do triângulo
    x = points[tri, 0]
    y = points[tri, 1]

    # Jacobiano e seu determinante
    J = np.array([[x[1] - x[0], x[2] - x[0]],
                  [y[1] - y[0], y[2] - y[0]]])
    detJ = np.linalg.det(J)
    invJ = np.linalg.inv(J)

    # Derivadas das funções de base no espaço real
    dphi_real = dphi_ref @ invJ

    # Matriz local S
    S_local = (dphi_real @ dphi_real.T) * detJ * w

    # Matriz local T
    T_local = np.full((3, 3), detJ * w / 12)
    np.fill_diagonal(T_local, detJ * w / 6)

    return S_local, T_local

# Função principal

def solve_gauss_1pt(radius=0.01, num_modes=6, filename='gauss1pt'):
    generate_circle_mesh(radius, filename)
    mesh = meshio.read(f"{filename}.msh")

    points = mesh.points[:, :2]
    triangles = mesh.cells_dict["triangle"]
    Nn = len(points)

    S = lil_matrix((Nn, Nn))
    T = lil_matrix((Nn, Nn))

    for tri in triangles:
        S_local, T_local = montar_matrizes_gauss_1pt(tri, points)
        for i in range(3):
            for j in range(3):
                S[tri[i], tri[j]] += S_local[i, j]
                T[tri[i], tri[j]] += T_local[i, j]

    vals, vecs = eigsh(S, k=num_modes, M=T, sigma=0, which='LM')

    kc = np.sqrt(np.real(vals))
    fc = 3e8 * kc / (2 * np.pi)

    return fc, vecs

# Execução

if __name__ == "__main__":
    fc, modos = solve_gauss_1pt()
    print("Frequências obtidas (GHz):", fc / 1e9)
