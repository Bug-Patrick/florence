import numpy as np
from Florence import *


def simple_laplace():
    """An example of solving the Laplace equation using
        fourth order hexahedral elements on a cube
    """

    # generate a linear hexahedral mesh on a cube
    mesh = Mesh()
    mesh.Cube(element_type="hex", nx=6, ny=6, nz=6)
    # generate the corresponding fourth order mesh
    mesh.GetHighOrderMesh(p=4)

    # set up boundary conditions
    def dirichlet_function(mesh):
        # create boundary flags - nan values would be treated as free boundary
        boundary_data = np.zeros(mesh.nnode)+np.nan
        # potential at left (Y=0)
        Y_0 = np.isclose(mesh.points[:,1],0)
        boundary_data[Y_0] = 0.
        # potential at right (Y=1)
        Y_1 = np.isclose(mesh.points[:,1],mesh.points[:,1].max())
        boundary_data[Y_1] = 10.

        return boundary_data

    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(dirichlet_function, mesh)

    # set up material
    material = IdealDielectric(mesh.InferSpatialDimension(), eps=2.35)
    # set up variational form
    formulation = LaplacianFormulation(mesh)
    # set up solver
    fem_solver = FEMSolver(optimise=True)
    # solve
    results = fem_solver.Solve( boundary_condition=boundary_condition,
                                material=material,
                                formulation=formulation,
                                mesh=mesh)

    # write results to vtk file
    # NB problem: /MeshGeneration/Mesh.py::GetLinearMesh(.) l.8214
    # inv[aranger].reshape(lmesh.nelem,nodeperelem) # aranger: index 216 is out of bounds for axis 0 with size 216
    # Solution attempt 1: aranger-1 is not enough -> maybe not number of elements = 216 passed to the function, but else?
    results.WriteVTK("laplacian_results")

    # plot results with florence
    # NB problem: /PostProcessing/PostProcess.py:: l.1662
    # sol = np.copy(self.sol[:mesh.nnode,:,:]) # IndexErro: too many indices for array: array is 2-dimensional, but 3 were indexed
    # maybe plot only works for 2D results? - but error says too many, when this 2D
    # results.Plot(configuration="deformed", quantity=0)



if __name__ == "__main__":
    simple_laplace()