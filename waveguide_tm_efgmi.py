
import subprocess

def generate_circle_mesh(radius=0.01, filename='circle'):
    geo_content = f"""
    lc = {radius/20};  // Tamanho característico da malha
    Point(1) = {{0, 0, 0, lc}};
    Point(2) = {{{radius}, 0, 0, lc}};
    Point(3) = {{0, {radius}, 0, lc}};
    Point(4) = {{-{radius}, 0, 0, lc}};
    Point(5) = {{0, -{radius}, 0, lc}};
    Circle(1) = {{2, 1, 3}};
    Circle(2) = {{3, 1, 4}};
    Circle(3) = {{4, 1, 5}};
    Circle(4) = {{5, 1, 2}};
    Line Loop(5) = {{1, 2, 3, 4}};
    Plane Surface(6) = {{5}};
    Physical Curve(100) = {{1, 2, 3, 4}};
    Physical Surface(200) = {{6}};
    Mesh 2;
    Save \"{filename}.msh\";
    """
    geo_file = f"{filename}.geo"
    with open(geo_file, "w") as f:
        f.write(geo_content)
    subprocess.run(["gmsh", geo_file, "-2", "-format", "msh2"], check=True)

import numpy as np
import meshio
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree

# Função peso de Schwarz singularizada
def weight_function(r, dmax, epsilon=1e-4):
    q = r / dmax
    w = ((1 - q)**4 * (4 * q + 1)) * (q < 1)
    return w * (q >= epsilon)

def construct_P(x):
    return np.array([1, x[0], x[1]])

def mls_shape_derivatives(xi, nodes, support_radius):
    tree = cKDTree(nodes)
    idx = tree.query_ball_point(xi, support_radius)
    neighbors = nodes[idx]
    w = weight_function(np.linalg.norm(neighbors - xi, axis=1), support_radius)

    if len(neighbors) < 3:
        return [], [], []

    P = np.array([construct_P(p) for p in neighbors])
    W = np.diag(w)
    A = P.T @ W @ P

    if np.linalg.cond(A) > 1e12:
        return [], [], []

    A_inv = np.linalg.inv(A)
    dN_dx, dN_dy = [], []
    p_x = np.array([0, 1, 0])
    p_y = np.array([0, 0, 1])
    for i, xj in enumerate(neighbors):
        pj = construct_P(xj)
        grad_phi = A_inv @ (P.T @ np.diag(w) @ (pj[:, None] - construct_P(xi)[:, None]))
        dN_dx.append(p_x @ grad_phi)
        dN_dy.append(p_y @ grad_phi)
    return idx, np.array(dN_dx), np.array(dN_dy)

def triangle_gauss_points(triangle):
    p1, p2, p3 = triangle
    area = 0.5 * abs(np.linalg.det([[p2[0] - p1[0], p2[1] - p1[1]],
                                    [p3[0] - p1[0], p3[1] - p1[1]]]))
    bary_coords = [(1/6, 1/6, 2/3), (1/6, 2/3, 1/6), (2/3, 1/6, 1/6)]
    gps = []
    weights = []
    for l1, l2, l3 in bary_coords:
        gp = l1 * p1 + l2 * p2 + l3 * p3
        gps.append(gp)
        weights.append(area / 3)
    return gps, weights

def solve_tm_efgmi(points, triangles, boundary_nodes, support_radius):
    Nn = len(points)
    node_id = np.ones(Nn, dtype=int)
    for i in boundary_nodes:
        node_id[i] = 0

    index = np.zeros(Nn, dtype=int)
    counter = 0
    for i in range(Nn):
        if node_id[i] == 1:
            counter += 1
            index[i] = counter

    Nf = counter
    S = lil_matrix((Nf, Nf))
    T = lil_matrix((Nf, Nf))

    for tri_idx in triangles:
        verts = points[tri_idx]
        gps, areas = triangle_gauss_points(verts)

        for gp, A in zip(gps, areas):
            idx, dphi_dx, dphi_dy = mls_shape_derivatives(gp, points, support_radius)
            if len(idx) == 0:
                continue
            for a in range(len(idx)):
                ia = idx[a]
                if node_id[ia] == 0: continue
                ia_idx = index[ia] - 1
                for b in range(len(idx)):
                    ib = idx[b]
                    if node_id[ib] == 0: continue
                    ib_idx = index[ib] - 1
                    Sa = (dphi_dx[a] * dphi_dx[b] + dphi_dy[a] * dphi_dy[b]) * A
                    Ta = A / len(idx)
                    S[ia_idx, ib_idx] += Sa
                    T[ia_idx, ib_idx] += Ta

    vals, vecs = eigsh(S, k=6, M=T, sigma=0, which='LM')
    vals = np.real(vals)
    modos = []
    for q in range(6):
        X0 = np.zeros(Nn)
        j = 0
        for i in range(Nn):
            if node_id[i] == 1:
                X0[i] = vecs[j, q]
                j += 1
        modos.append(X0)
    return np.sqrt(vals), modos

# MAIN

if __name__ == "__main__":
    generate_circle_mesh(radius=0.01, filename="circle")
    mesh = meshio.read("circle.msh")
    points = mesh.points[:, :2]
    triangles = mesh.cells_dict["triangle"]

    mesh = meshio.read("circle.msh")
    points = mesh.points[:, :2]
    triangles = mesh.cells_dict["triangle"]

    boundary_nodes = set()
    for cell_block, phys_ids in zip(mesh.cells, mesh.cell_data_dict["gmsh:physical"].values()):
        if cell_block.type == "line":
            for i, line in enumerate(cell_block.data):
                if phys_ids[i] == 100:
                    boundary_nodes.update(line)

    kc_vals, modos = solve_tm_efgmi(points, triangles, boundary_nodes, support_radius=0.02)
    print("Valores de kc:")
    print(kc_vals)

    # Plot do primeiro modo
    plt.figure(figsize=(6,5))
    plt.tricontourf(points[:,0], points[:,1], triangles, modos[0], levels=100, cmap='jet')
    plt.colorbar(label='Amplitude')
    plt.title("Modo 1 - TM via EFGMI")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig("modo1_tm_efgmi.png", dpi=300)
    plt.show()
