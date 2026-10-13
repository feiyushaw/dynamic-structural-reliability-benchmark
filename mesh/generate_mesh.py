"""
Mesh generator for Cantilever Steel Beam.
Generates structured Hex8 mesh consistent with geometry/cantilever_beam.geo:
Dimensions: L = 300 mm, H = 40 mm, B = 10 mm.
Mesh divisions: 60 x 8 x 2 (960 elements, 5mm resolution).
"""

import numpy as np
from mpi4py import MPI
import dolfinx
from dolfinx import mesh


def create_cantilever_mesh(
    comm=MPI.COMM_WORLD,
    length: float = 0.3,
    height: float = 0.04,
    thickness: float = 0.01,
    nx: int = 60,
    ny: int = 8,
    nz: int = 2,
):
    """Create 3D structured Hex8 mesh using DOLFINx."""
    domain = mesh.create_box(
        comm,
        [np.array([0.0, 0.0, 0.0]), np.array([length, height, thickness])],
        [nx, ny, nz],
        cell_type=dolfinx.mesh.CellType.hexahedron,
    )
    return domain


if __name__ == "__main__":
    domain = create_cantilever_mesh()
    num_cells = domain.topology.index_map(domain.topology.dim).size_local
    print(f"Cantilever beam mesh generated successfully: {num_cells} Hex8 elements.")
