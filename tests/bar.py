import numpy as np
from Florence import *


def bar_problem_setup(increments=1, force_direction=1, force_magnitude=-100):
    """Setup bar problem common between VariationalFormulations and Materials
        increments: number of load increments for nonlinear analysis
        force_direction: 0 for x, 1 for y, 2 for z direction; standard is y direction!
        force_magnitude: magnitude of the applied force in the specified direction
    """

    # read a mesh from a gmsh file
    mesh = Mesh()
    mesh.ReadGmsh(os.path.join(PWD(__file__), "bar.msh"), element_type="tet")
    mesh.ndim = mesh.InferSpatialDimension()

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

        # at right (X=0, Y=Z=0.5) a x-direction force of -100
        NeumannBC = (0., 0.5, 0.5)
        idx = np.where((mesh.points[:] == NeumannBC).all(axis=1))[0]
        boundary_data[idx, force_direction] = force_magnitude
        return boundary_data

    increment_step = increments
    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(DirichletFunc, mesh)
    boundary_condition.SetNeumannCriteria(NeumannFuncDyn, mesh)

    # set up solver
    # careful to not enable low-level dispatcher: if(has_low_level_dispatcher != optimise): has_low_level_dispatcher = True
    fem_solver = FEMSolver(
        number_of_load_increments=increment_step,
        analysis_type="static",
        # analysis_subtype="explicit", # Explicit or implicit??
        analysis_nature="nonlinear",
        # optimise=True, # has_low_level_dispatcher=False, # True-False is bad combination: RuntimeError: Cannot dispatch to low level module since material NeoHookeanF does not support it
        memory_store_frequency=20)
    
    return mesh, boundary_condition, fem_solver





def bar_NL_tests():
    """An use case of solving a bar problem using
        linear elements read from a gmsh file

        Working file, testing framework
    """

    # Read gmsh file, create boundary conditions and solver
    mesh, boundary_condition, fem_solver = bar_problem_setup()

    # Set material data
    youngs_modulus = 502000
    poissons_ratio = 0.4
    # material = NeoHookean(mesh.ndim, youngs_modulus=502000, poissons_ratio=0.4)
    lame_parameter_1 = youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))
    lame_parameter_2 = youngs_modulus / (2 * (1 + poissons_ratio))

    material = NeoHookeanF(mesh.ndim, lamb=lame_parameter_1, mu=lame_parameter_2)  # first parameter ndim
    # NeoHookeanF ? OgdenNeoHookeanC ? StVenantKirchhoffC


    # set up variational form
    formulation = FBasedDisplacementFormulation(mesh)
    # FBasedDisplacmentFormulation; CBasedDisplacementFormulation; DisplacementFormulation

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



def bar_MR(simulation_type="F", stabilise_tangents=True):
    """An use case of solving a bar problem using
        linear elements read from a gmsh file

        simulation_type: F or TL for FBased or standard Total Lagrangian formulation
        stabilise_tangents: whether to stabilise tangents 
    """

    # Read gmsh file, create boundary conditions and solver
    mesh, boundary_condition, fem_solver = bar_problem_setup(increments=1, force_direction=1, force_magnitude=-100000)

    # Set material data
    youngs_modulus = 502000
    poissons_ratio = 0.4
    # material = NeoHookean(mesh.ndim, youngs_modulus=502000, poissons_ratio=0.4)
    lamb = youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))
    mu = youngs_modulus / (2 * (1 + poissons_ratio))
    # lamb, mu = 717142.8571428574, 179285.7142857143

    # split mu1=C10 and mu2=C01 for Mooney-Rivlin?
    # Gemini says: C10 = 0.41 MPa and C01 = 0.43 MPa? were to gain reliable values?

    if(simulation_type == "F"):
        # Set material data
        material = MooneyRivlinF(mesh.ndim, lamb=lamb, mu1=mu, mu2=mu, minJ=1, stabilise_tangents=stabilise_tangents)

        # set up variational form
        formulation = FBasedDisplacementFormulation(mesh)
    else:
        # Set material data
        material = MooneyRivlin(mesh.ndim, lamb=lamb, mu1=mu, mu2=mu, minJ=1)

        # set up variational form
        formulation = DisplacementFormulation(mesh)

    solution = fem_solver.Solve(formulation=formulation, material=material, mesh=mesh,
        boundary_condition=boundary_condition)

    # check validity ?
    solution_vectors = solution.GetSolutionVectors()

    # export 0.result field to vtk file
    solution.WriteVTK("bar_MR_"+simulation_type, quantity=0)




def bar_NH_F():
    """An use case of solving a bar problem using
        linear elements read from a gmsh file

        Inspired by Car Crash analysis example but removing contact and dynamics
    """

    # Read gmsh file, create boundary conditions and solver
    mesh, boundary_condition, fem_solver = bar_problem_setup()

    # Set material data
    material = OgdenNeoHookeanF(mesh.ndim, youngs_modulus=502000, poissons_ratio=0.4, minJ=1)

    # set up variational form
    formulation = FBasedDisplacementFormulation(mesh)

    solution = fem_solver.Solve(formulation=formulation, material=material, mesh=mesh,
        boundary_condition=boundary_condition)

    # check validity ?
    solution_vectors = solution.GetSolutionVectors()

    # export 0.result field to vtk file
    solution.WriteVTK("bar_NH_F", quantity=0)





def bar_NH():
    """An use case of solving a bar problem using
        linear elements read from a gmsh file
    """

    # Read gmsh file, create boundary conditions and solver
    mesh, boundary_condition, fem_solver = bar_problem_setup()

    # Set material data
    material = NeoHookean(mesh.ndim, youngs_modulus=502000, poissons_ratio=0.4)

    # set up variational form
    formulation = DisplacementFormulation(mesh)

    solution = fem_solver.Solve(formulation=formulation, material=material, mesh=mesh,
        boundary_condition=boundary_condition)

    # check validity ?
    solution_vectors = solution.GetSolutionVectors()

    # export 1.result field (0-2:u_x,u_y,u_z) to vtk file
    solution.WriteVTK("bar_NH", quantity=1)

    return solution_vectors






if __name__ == "__main__":
    bar_MR(simulation_type="F", stabilise_tangents=True)
    bar_MR(simulation_type="TL", stabilise_tangents=False)
    #bar_NH()
    #bar_NH_F()
    #bar_NL_tests()