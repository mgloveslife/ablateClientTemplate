static char help[] = "Example ablate Client using incompressible flow and native PETSc api";

/** Example Arguments
 -dm_plex_separate_marker -dm_refine 0
 -vel_petscspace_degree 2 -pres_petscspace_degree 1 -temp_petscspace_degree 1
 -dmts_check .001 -ts_max_steps 4 -ts_dt 0.1
 -ksp_type fgmres -ksp_gmres_restart 10 -ksp_rtol 1.0e-9 -ksp_error_if_not_converged
 -pc_type fieldsplit -pc_fieldsplit_0_fields 0,2 -pc_fieldsplit_1_fields 1 -pc_fieldsplit_type schur -pc_fieldsplit_schur_factorization_type full
 -fieldsplit_0_pc_type lu
 -fieldsplit_pressure_ksp_rtol 1e-10 -fieldsplit_pressure_pc_type jacobi
 */

#include <petsc.h>
#include <environment/runEnvironment.hpp>
#include <monitors/hdf5Monitor.hpp>
#include <parameters/mapParameters.hpp>
#include "flow/boundaryConditions/essential.hpp"
#include "flow/incompressibleFlow.hpp"
#include "mathFunctions/functionFactory.hpp"
#include "mesh/dmWrapper.hpp"
#include "parameters/petscOptionParameters.hpp"
#include "utilities/petscError.hpp"

using namespace ablate;

typedef PetscErrorCode (*ExactFunction)(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);

typedef void (*IntegrandTestFunction)(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt *uOff, const PetscInt *uOff_x, const PetscScalar *u, const PetscScalar *u_t, const PetscScalar *u_x,
                                      const PetscInt *aOff, const PetscInt *aOff_x, const PetscScalar *a, const PetscScalar *a_t, const PetscScalar *a_x, PetscReal t, const PetscReal *X,
                                      PetscInt numConstants, const PetscScalar *constants, PetscScalar *f0);

#define SourceFunction(FUNC)            \
    FUNC(PetscInt dim,                  \
         PetscInt Nf,                   \
         PetscInt NfAux,                \
         const PetscInt uOff[],         \
         const PetscInt uOff_x[],       \
         const PetscScalar u[],         \
         const PetscScalar u_t[],       \
         const PetscScalar u_x[],       \
         const PetscInt aOff[],         \
         const PetscInt aOff_x[],       \
         const PetscScalar a[],         \
         const PetscScalar a_t[],       \
         const PetscScalar a_x[],       \
         PetscReal t,                   \
         const PetscReal X[],           \
         PetscInt numConstants,         \
         const PetscScalar constants[], \
         PetscScalar f0[])

// store the pointer to the provided test function from the solver
static IntegrandTestFunction f0_v_original;
static IntegrandTestFunction f0_w_original;
static IntegrandTestFunction f0_q_original;

static PetscErrorCode SetInitialConditions(TS ts, Vec u) {
    DM dm;
    PetscReal t;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = TSGetDM(ts, &dm);
    CHKERRQ(ierr);
    ierr = TSGetTime(ts, &t);
    CHKERRQ(ierr);

    // This function Tags the u vector as the exact solution.  We need to copy the values to prevent this.
    Vec e;
    ierr = VecDuplicate(u, &e);
    CHKERRQ(ierr);
    ierr = DMComputeExactSolution(dm, t, e, NULL);
    CHKERRQ(ierr);
    ierr = VecCopy(e, u);
    CHKERRQ(ierr);
    ierr = VecDestroy(&e);
    CHKERRQ(ierr);

    // Get the flowData
    ablate::flow::Flow *flow;
    ierr = DMGetApplicationContext(dm, &flow);
    CHKERRQ(ierr);

    // get the flow to apply the completeFlowInitialization method
    flow->CompleteFlowInitialization(dm, u);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode MonitorErrorCustom(TS ts, PetscInt step, PetscReal crtime, Vec u, void *ctx) {
    PetscErrorCode (*exactFuncs[3])(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);
    void *ctxs[3];
    DM dm;
    PetscDS ds;
    Vec v;
    PetscReal ferrors[3];
    PetscInt f;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = TSGetDM(ts, &dm);
    CHKERRQ(ierr);
    ierr = DMGetDS(dm, &ds);
    CHKERRQ(ierr);

    for (f = 0; f < 3; ++f) {
        ierr = PetscDSGetExactSolution(ds, f, &exactFuncs[f], &ctxs[f]);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
    }
    ierr = DMComputeL2FieldDiff(dm, crtime, exactFuncs, ctxs, u, ferrors);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Timestep: %04d time = %-8.4g \t L_2 Error: [%2.3g, %2.3g, %2.3g]\n", (int)step, (double)crtime, (double)ferrors[0], (double)ferrors[1], (double)ferrors[2]);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    PetscFunctionReturn(0);
}

// helper functions for generated code
static PetscReal Power(PetscReal x, PetscInt exp) { return PetscPowReal(x, exp); }
static PetscReal Cos(PetscReal x) { return PetscCosReal(x); }
static PetscReal Sin(PetscReal x) { return PetscSinReal(x); }

/*
  CASE: incompressible quadratic
  In 2D we use exact solution:

    u = t + x^2 + y^2
    v = t + 2x^2 - 2xy
    p = x + y - 1
    T = t + x + y
  so that

    \nabla \cdot u = 2x - 2x = 0

  see docs/content/formulations/incompressibleFlow/solutions/Incompressible_2D_Quadratic_MMS.nb
*/
static PetscErrorCode incompressible_quadratic_u(PetscInt Dim, PetscReal time, const PetscReal *X, PetscInt Nf, PetscScalar *u, void *ctx) {
    u[0] = time + X[0] * X[0] + X[1] * X[1];
    u[1] = time + 2.0 * X[0] * X[0] - 2.0 * X[0] * X[1];
    return 0;
}
static PetscErrorCode incompressible_quadratic_u_t(PetscInt Dim, PetscReal time, const PetscReal *X, PetscInt Nf, PetscScalar *u, void *ctx) {
    u[0] = 1.0;
    u[1] = 1.0;
    return 0;
}

static PetscErrorCode incompressible_quadratic_p(PetscInt Dim, PetscReal time, const PetscReal *X, PetscInt Nf, PetscScalar *p, void *ctx) {
    p[0] = X[0] + X[1] - 1.0;
    return 0;
}

static PetscErrorCode incompressible_quadratic_T(PetscInt Dim, PetscReal time, const PetscReal *X, PetscInt Nf, PetscScalar *T, void *ctx) {
    T[0] = time + X[0] + X[1];
    return 0;
}
static PetscErrorCode incompressible_quadratic_T_t(PetscInt Dim, PetscReal time, const PetscReal *X, PetscInt Nf, PetscScalar *T, void *ctx) {
    T[0] = 1.0;
    return 0;
}

/* f0_v = du/dt - f */
static void SourceFunction(f0_incompressible_quadratic_v) {
    f0_v_original(dim, Nf, NfAux, uOff, uOff_x, u, u_t, u_x, aOff, aOff_x, a, a_t, a_x, t, X, numConstants, constants, f0);
    const PetscReal rho = 1.0;
    const PetscReal S = constants[0];   // STROUHAL
    const PetscReal mu = constants[3];  // MU
    const PetscReal R = constants[1];   // REYNOLDS
    const PetscReal x = X[0];
    const PetscReal y = X[1];

    f0[0] -= 1 - (4. * mu) / R + rho * S + 2 * rho * y * (t + 2 * Power(x, 2) - 2 * x * y) + 2 * rho * x * (t + Power(x, 2) + Power(y, 2));
    f0[1] -= 1 - (4. * mu) / R + rho * S - 2 * rho * x * (t + 2 * Power(x, 2) - 2 * x * y) + rho * (4 * x - 2 * y) * (t + Power(x, 2) + Power(y, 2));
}

/* f0_w = dT/dt + u.grad(T) - Q */
static void SourceFunction(f0_incompressible_quadratic_w) {
    f0_w_original(dim, Nf, NfAux, uOff, uOff_x, u, u_t, u_x, aOff, aOff_x, a, a_t, a_x, t, X, numConstants, constants, f0);

    const PetscReal rho = 1.0;
    const PetscReal S = constants[0];   // STROUHAL
    const PetscReal Cp = constants[5];  // CP
    const PetscReal x = X[0];
    const PetscReal y = X[1];

    f0[0] -= Cp * rho * (S + 2 * t + 3 * Power(x, 2) - 2 * x * y + Power(y, 2));
}

int main(int argc, char *argv[]) {
    DM dm; /* problem definition */
    TS ts; /* timestepper */
    PetscReal t;
    PetscErrorCode ierr;

    // initialize petsc and mpi
    PetscInitialize(&argc, &argv, NULL, help);
    {
        // setup the run environment
        ablate::parameters::MapParameters runEnvironmentParameters(std::map<std::string, std::string>{{"title", "clientExample"}});
        ablate::environment::RunEnvironment::Setup(runEnvironmentParameters);

        // setup the ts
        TSCreate(PETSC_COMM_WORLD, &ts) >> ablate::checkError;
        // Create a mesh
        // hard code the problem setup
        DMPlexCreateBoxMesh(PETSC_COMM_WORLD, 2, PETSC_FALSE, NULL, NULL, NULL, NULL, PETSC_TRUE, &dm) >> ablate::checkError;
        DMSetFromOptions(dm) >> ablate::checkError;
        TSSetDM(ts, dm) >> ablate::checkError;
        TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP) >> ablate::checkError;

        // pull the parameters from the petsc options
        auto parameters = std::make_shared<ablate::parameters::PetscOptionParameters>();

        // set the init and exact solution from c functions
        auto velocityExact = std::make_shared<mathFunctions::FieldSolution>("velocity", mathFunctions::Create(incompressible_quadratic_u), mathFunctions::Create(incompressible_quadratic_u_t));
        auto pressureExact = std::make_shared<mathFunctions::FieldSolution>("pressure", mathFunctions::Create(incompressible_quadratic_p));
        auto temperatureExact = std::make_shared<mathFunctions::FieldSolution>("temperature", mathFunctions::Create(incompressible_quadratic_T), mathFunctions::Create(incompressible_quadratic_T_t));

        // Create the flow object
        auto flowObject = std::make_shared<flow::IncompressibleFlow>(
            "testFlow",
            std::make_shared<mesh::DMWrapper>(dm),
            parameters,
            nullptr,
            /* initialization functions */
            std::vector<std::shared_ptr<mathFunctions::FieldSolution>>{velocityExact, pressureExact, temperatureExact},
            /* boundary conditions */
            std::vector<std::shared_ptr<flow::boundaryConditions::BoundaryCondition>>{
                std::make_shared<flow::boundaryConditions::Essential>(
                    "velocity", "velocity wall", "marker", std::vector<int>{3, 1, 2, 4}, mathFunctions::Create(incompressible_quadratic_u), mathFunctions::Create(incompressible_quadratic_u_t)),
                std::make_shared<flow::boundaryConditions::Essential>(
                    "temperature", "temp wall", "marker", std::vector<int>{3, 1, 2, 4}, mathFunctions::Create(incompressible_quadratic_T), mathFunctions::Create(incompressible_quadratic_T_t)),
            },
            /* aux field updates */
            std::vector<std::shared_ptr<mathFunctions::FieldSolution>>{},
            /* exact */
            std::vector<std::shared_ptr<mathFunctions::FieldSolution>>{velocityExact, pressureExact, temperatureExact});

        // Override problem with source terms, boundary, and set the exact solution
        {
            PetscDS prob;
            DMGetDS(dm, &prob) >> ablate::checkError;

            // V, W Test Function
            IntegrandTestFunction tempFunctionPointer;
            PetscDSGetResidual(prob, 0, &f0_v_original, &tempFunctionPointer) >> ablate::checkError;
            PetscDSSetResidual(prob, 0, f0_incompressible_quadratic_v, tempFunctionPointer) >> ablate::checkError;

            PetscDSGetResidual(prob, 2, &f0_w_original, &tempFunctionPointer) >> ablate::checkError;
            PetscDSSetResidual(prob, 2, f0_incompressible_quadratic_w, tempFunctionPointer) >> ablate::checkError;
        }
        flowObject->CompleteProblemSetup(ts);

        // Setup the TS
        TSSetFromOptions(ts) >> ablate::checkError;

        // Set initial conditions from the exact solution
        TSSetComputeInitialCondition(ts, SetInitialConditions) >> ablate::checkError;

        // use a monitor from the ablate lib
        auto hdf5Monitor = std::make_shared<ablate::monitors::Hdf5Monitor>();
        hdf5Monitor->Register(flowObject);
        TSMonitorSet(ts, hdf5Monitor->GetPetscFunction(), hdf5Monitor->GetContext(), NULL) >> ablate::checkError;

        // Add in a custom monitor
        TSMonitorSet(ts, MonitorErrorCustom, NULL, NULL) >> ablate::checkError;

        TSSolve(ts, flowObject->GetSolutionVector()) >> ablate::checkError;

        // Compare the actual vs expected values
        DMTSCheckFromOptions(ts, flowObject->GetSolutionVector()) >> ablate::checkError;

        // Cleanup
        DMDestroy(&dm) >> ablate::checkError;
        TSDestroy(&ts) >> ablate::checkError;
    }
    return PetscFinalize();
}