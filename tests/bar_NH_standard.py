import numpy as np
from Florence import *

from bar import bar_NH


def bar():
    """An use case of solving a bar problem using
        linear elements read from a gmsh file

        Inspired by Car Crash analysis example but removing contact and dynamics
    """

    # read a mesh from a gmsh file
    mesh = Mesh()
    mesh.ReadGmsh(os.path.join(PWD(__file__),"bar.msh"),element_type="tet")
    mesh.ndim = mesh.InferSpatialDimension()

    # Set material data
    material = NeoHookean(mesh.ndim, youngs_modulus=502000, poissons_ratio=0.4)

    def DirichletFunc(mesh):
        # homogenous Dirichlet boundary at elements 10 11 22 23 49 - nan values as free boundary
        boundary_data = np.zeros((mesh.nnode,3))+np.nan
        # at left (X=5)
        X_0 = np.isclose(mesh.points[:,0],5.)
        boundary_data[X_0,:] = (0., 0., 0.)
        return boundary_data


    def NeumannFuncDyn(mesh):
        # Neumann boundary with -100 force in x or y direction at element 44 - nan values as free boundary
        boundary_data = np.zeros((mesh.points.shape[0],3))+np.nan
        
        # at right (X=0, Y=Z=0.5) a x-direction force of -100
        NeumannBC = (0., 0.5, 0.5)
        idx = np.where((mesh.points[:] == NeumannBC).all(axis=1))[0]
        boundary_data[idx,0] = -100.
        return boundary_data


    increment_step = 1
    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(DirichletFunc, mesh)
    boundary_condition.SetNeumannCriteria(NeumannFuncDyn, mesh)

    # set up variational form
    formulation = DisplacementFormulation(mesh)

    # set up solver
    fem_solver = FEMSolver(
        number_of_load_increments=increment_step,
        analysis_type="static",
        # analysis_subtype="explicit", # Explicit or implicit??
        analysis_nature="nonlinear",
        optimise=True,
        memory_store_frequency=20)

    solution = fem_solver.Solve(formulation=formulation, material=material, mesh=mesh,
        boundary_condition=boundary_condition)

    # check validity ?
    solution_vectors = solution.GetSolutionVectors()

    # export 0.result field to vtk file
    solution.WriteVTK("bar_NH", quantity=0)



if __name__ == "__main__":
    # bar()
    bar_NH()