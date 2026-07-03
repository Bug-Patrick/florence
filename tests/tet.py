import numpy as np
from Florence import *


def tet_problem_setup(stabilise=False,increments=1, force_direction=1, force_magnitude=0.001):
    """Setup bar problem common between VariationalFormulations and Materials
        increments: number of load increments for nonlinear analysis
        force_direction: 0 for x, 1 for y, 2 for z direction; standard is y direction!
        force_magnitude: magnitude of the applied force in the specified direction
    """

    # read a mesh from a gmsh file
    mesh = Mesh()
    mesh.ReadGmsh(os.path.join(PWD(__file__), "tet.gmsh"), element_type="tet")
    mesh.ndim = mesh.InferSpatialDimension()


    def DirichletFunc(mesh):
        # homogenous Dirichlet boundary at nodes 0 1 2 - nan values as free boundary
        boundary_data = np.zeros((mesh.nnode, 3)) + np.nan
        print(mesh.points)
        # at nodes 0:2 al DoF : are fixed
        boundary_data[0:2, :] = (0., 0., 0.)
        return boundary_data

    def NeumannFuncDyn(mesh):
        # Neumann boundary with 0.001 force in zdirection at node 3 - nan values as free boundary
        boundary_data = np.zeros((mesh.points.shape[0], 3)) + np.nan

        # at node 3 a z(x_3)-direction force of 0.001
        boundary_data[3, 2] = 0.001
        return boundary_data

    increment_step = increments
    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(DirichletFunc, mesh)
    boundary_condition.SetNeumannCriteria(NeumannFuncDyn, mesh)

    # set up solver
    # careful to not enable low-level dispatcher: if(has_low_level_dispatcher != optimise): has_low_level_dispatcher = True
    fem_solver = FEMSolver(
        # === Basic setup ===
        number_of_load_increments=increment_step,
        analysis_type="static",
        # analysis_subtype="explicit", # Explicit or implicit?? for dynamics?
        analysis_nature="nonlinear",
        # === Research code setup ===
        stabilise_local_system=stabilise, # stabilise using analytic eigensystem
        #stabilise_global_system # stabilise using eigenvalue decomposition?
        # === Code optimisation ===
        # optimise=True, # has_low_level_dispatcher=False, # True-False is bad combination: RuntimeError: Cannot dispatch to low level module since material NeoHookeanF does not support it
        # === Debugging prints ===
        print_incremental_log=True,
        save_incremental_solution=True,
        incremental_solution_filename="florence Tet Sol Incr",
        incremental_solution_save_frequency=1,
        break_at_increment=0,
        memory_store_frequency=1#,
        # === Nonlinear solver advanced settings ===
        #activate_line_search
        #activate_arc_length
        #newton_raphson_tolerance
        #newton_raphson_solution_tolerance
        #maximum_iteration_for_newton_raphson
        #nonlinear_iterative_technique
        #line_search_technique
        )
    
    return mesh, boundary_condition, fem_solver





def tet_MR(simulation_type="F", stabilise_tangents=True):
    """An use case of solving a bar problem using
        linear elements read from a gmsh file

        simulation_type: F or TL for FBased or standard Total Lagrangian formulation
        stabilise_tangents: whether to stabilise tangents 
    """

    # Read gmsh file, create boundary conditions and solver
    mesh, boundary_condition, fem_solver = tet_problem_setup(stabilise=stabilise_tangents,increments=1, force_direction=1, force_magnitude=100)

    # Set material data
    youngs_modulus = 502000
    poissons_ratio = 0.4
    lamb = youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))
    mu = youngs_modulus / (2 * (1 + poissons_ratio))
    # lamb, mu = 717142.8571428574, 179285.7142857143

    # split mu1=C10 and mu2=C01 for Mooney-Rivlin?
    # Gemini says: C10 = 0.41 MPa and C01 = 0.43 MPa? are reliable values?

    if(simulation_type == "F"):
        # Set material data
        material = MooneyRivlinF(mesh.ndim, lamb=lamb, mu1=mu, mu2=mu, minJ=0.5, stabilise_tangents=stabilise_tangents)

        # set up variational form
        print("Stabilisation: "+str(stabilise_tangents))
        formulation = FBasedDisplacementFormulation(mesh)
    else:
        # Set material data
        material = MooneyRivlin(mesh.ndim, lamb=lamb, mu1=mu, mu2=mu, minJ=0.5)

        # set up variational form
        formulation = DisplacementFormulation(mesh)

    solution = fem_solver.Solve(formulation=formulation, material=material, mesh=mesh,
        boundary_condition=boundary_condition)

    # check validity ?
    solution_vectors = solution.GetSolutionVectors()

    # export 0.result field to vtk file
    solution.WriteVTK("tet_MR_" + simulation_type, quantity=0)
    solution.WriteVTK("tet_MR_" + simulation_type, quantity=1)
    solution.WriteVTK("tet_MR_" + simulation_type, quantity=2)






if __name__ == "__main__":
    tet_MR(simulation_type="F", stabilise_tangents=False) # Validate unstabilised stiffness!
    # bar_MR(simulation_type="F", stabilise_tangents=True)
    # bar_MR(simulation_type="TL", stabilise_tangents=False)
    # bar_NH(simulation_type="F", material_formulation=3, stabilise_tangents=True)
    # bar_NH(simulation_type="TL", stabilise_tangents=False)
    #bar_NL_tests()