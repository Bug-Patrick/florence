import numpy as np
from Florence import *


def bar_NL_tests():
    """An use case of solving a bar problem using
        linear elements read from a gmsh file

        Working file, testing framework
        Inspired by Car Crash analysis example but removing contact and dynamics
    """

    # read a mesh from a gmsh file
    mesh = Mesh()
    mesh.ReadGmsh(os.path.join(PWD(__file__), "bar.msh"), element_type="tet")
    mesh.ndim = mesh.InferSpatialDimension()

    # Set material data
    youngs_modulus = 502000
    poissons_ratio = 0.4
    # material = NeoHookean(mesh.ndim, youngs_modulus=502000, poissons_ratio=0.4)
    lame_parameter_1 = youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))
    lame_parameter_2 = youngs_modulus / (2 * (1 + poissons_ratio))

    material = NeoHookean(mesh.ndim, lamb=lame_parameter_1, mu=lame_parameter_2)  # first parameter ndim
    # NeoHookeanF ? OgdenNeoHookeanC ? StVenantKirchhoffC

    def DirichletFunc(mesh):
        # homogenous Dirichlet boundary at elements 10 11 22 23 49 - nan values as free boundary
        boundary_data = np.zeros((mesh.nnode, 3)) + np.nan
        # at left (X=5)
        X_0 = np.isclose(mesh.points[:, 0], 5.)
        boundary_data[X_0, :] = (0., 0., 0.)
        return boundary_data

    def NeumannFuncDyn(mesh):
        # Neumann boundary with -100 force in x or y direction at element 44 - nan values as free boundary
        boundary_data = np.zeros((mesh.points.shape[0], 3)) + np.nan

        # at right (X=0, Y=Z=0.5)
        NeumannBC = (0., 0.5, 0.5)
        idx = np.where((mesh.points[:] == NeumannBC).all(axis=1))[0]
        boundary_data[idx, 0] = -100.
        # boundary_data[idx,1] = -100.
        # boundary_data[idx,2] = 0.
        return boundary_data

    increment_step = 1
    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(DirichletFunc, mesh)
    boundary_condition.SetNeumannCriteria(NeumannFuncDyn, mesh)

    # provide QuadratureRule and FunctionSpace for VariationalFormulation?
    

    # set up variational form
    formulation = CBasedDisplacementFormulation(mesh)
    # FBasedDisplacmentFormulation; CBasedDisplacementFormulation; DisplacementFormulation

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

    # write results to vtk file
    solution.WriteVTK("bar_work", quantity=0)
    # Function from PostProcess.py:def WriteVTK(.): quantity is index for the field to be plotted; 'all' for every field
    # NB problem for laplacian example: /MeshGeneration/Mesh.py::GetLinearMesh(.) l.8214
    # Solution: aranger should be a range of all elements not all nodes of the elements ~> is this fix for all?

    # plot results with florence
    # NB problem: /PostProcessing/PostProcess.py:: l.1662
    # sol = np.copy(self.sol[:mesh.nnode,:,:]) # IndexErro: too many indices for array: array is 2-dimensional, but 3 were indexed
    # maybe plot only works for 2D results? - but error says too many, when this 2D
    # results.Plot(configuration="deformed", quantity=0)


def bar_NH():
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
    bar_NH()
    bar_NL_tests()