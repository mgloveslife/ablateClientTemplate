#include </Users/jjmarzia/Desktop/ablate/petsc/include/petsc.h>
#include </Users/jjmarzia/Desktop/ablate/ablate/src/environment/runEnvironment.hpp>
#include </Users/jjmarzia/Desktop/ablate/ablate/src/mathFunctions/functionFactory.hpp>
#include <memory>
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/builder.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/domain/boxMesh.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/domain/modifiers/distributeWithGhostCells.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/domain/modifiers/ghostBoundaryCells.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/eos/perfectGas.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/finiteVolume/boundaryConditions/ghost.hpp"

#include "/Users/jjmarzia/Desktop/ablate/ablate/src/finiteVolume/compressibleFlowSolver.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/finiteVolume/fluxCalculator/ausm.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/io/interval/fixedInterval.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/monitors/curveMonitor.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/monitors/timeStepMonitor.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/parameters/mapParameters.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/utilities/petscUtilities.hpp"

#include "/Users/jjmarzia/Desktop/ablate/ablate/src/io/hdf5MultiFileSerializer.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/io/interval/simulationTimeInterval.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/finiteVolume/compressibleFlowFields.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/finiteVolume/fieldFunctions/euler.hpp"
#include "/Users/jjmarzia/Desktop/ablate/ablate/src/finiteVolume/fieldFunctions/compressibleFlowState.hpp"

typedef struct {
    PetscReal gamma;
    PetscReal Rgas;
    PetscReal temperature;
    PetscReal pressure;
    PetscReal volumeFraction;
    PetscReal rho;


//    PetscReal length;
//    PetscReal rhoL;
//    PetscReal uL;
//    PetscReal pL;
//    PetscReal rhoR;
//    PetscReal uR;
//    PetscReal pR;
} InitialConditions;

static PetscErrorCode SetInitialCondition(PetscInt dim, PetscReal time, const PetscReal x[], const PetscReal y[], PetscInt Nf, PetscScalar *u, void *ctx) {
    InitialConditions *initialConditions = (InitialConditions *)ctx;

//    if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1) {
//        initialConditions->pressure = 1e5 + 0.251036495228486;
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.1025) {
//        initialConditions->pressure = 1e5 + 0.209197079357072;
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.21) {
//        initialConditions->pressure = 1e5+0.167357663485657;           // initialize pressure
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.3225) {
//        initialConditions->pressure = 1e5+0.125518247614243;
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.44) {
//        initialConditions->pressure = 1e5 + 0.083678831742829;
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.5625) {
//        initialConditions->pressure = 1e5+0.041839415871414;
//    } else {
//        initialConditions->pressure = 1e5;
//    }
//
//    if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1) {
//        initialConditions->volumeFraction = 0;
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.1025) {
//        initialConditions->volumeFraction = 0.166666666666667;
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.21) {
//        initialConditions->volumeFraction = 0.333333333333333;           // initialize vol fraction
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.3225) {
//        initialConditions->volumeFraction = 0.5;
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.44) {
//        initialConditions->volumeFraction = 0.666666666666667;
//    } else if (PetscSqr(x[0]) + PetscSqr(y[0]) < 1.5625) {
//        initialConditions->volumeFraction = 0.833333333333333;
//    } else {
//        initialConditions->volumeFraction = 1.0;
//    }

    PetscReal u = 0;
    PetscReal v = 0;
    PetscReal e = initialConditions->pressure / ((initialConditions->gamma - 1.0) * initialConditions->rho);
    PetscReal et = e + 0.5 * u * u;

    u[ablate::finiteVolume::CompressibleFlowFields::RHOU] = initialConditions->rho*u;
    u[ablate::finiteVolume::CompressibleFlowFields::RHOV] = initialConditions->rho*v;

//    if (x[0] < initialConditions->length / 2.0) {
//        u[ablate::finiteVolume::CompressibleFlowFields::RHO] = initialConditions->rhoL;
//        u[ablate::finiteVolume::CompressibleFlowFields::RHOU + 0] = initialConditions->rhoL * initialConditions->uL;
//
//        PetscReal e = initialConditions->pL / ((initialConditions->gamma - 1.0) * initialConditions->rhoL);
//        PetscReal et = e + 0.5 * PetscSqr(initialConditions->uL);
//        u[ablate::finiteVolume::CompressibleFlowFields::RHOE] = et * initialConditions->rhoL;
//
//    } else {
//        u[ablate::finiteVolume::CompressibleFlowFields::RHO] = initialConditions->rhoR;
//        u[ablate::finiteVolume::CompressibleFlowFields::RHOU + 0] = initialConditions->rhoR * initialConditions->uR;
//
//        PetscReal e = initialConditions->pR / ((initialConditions->gamma - 1.0) * initialConditions->rhoR);
//        PetscReal et = e + 0.5 * PetscSqr(initialConditions->uR);
//        u[ablate::finiteVolume::CompressibleFlowFields::RHOE] = et * initialConditions->rhoR;
//    }

    return 0;
}

static PetscErrorCode PhysicsBoundary_Euler(PetscReal time, const PetscReal *c, const PetscReal *n, const PetscScalar *a_xI, PetscScalar *a_xG, void *ctx) {
    InitialConditions *initialConditions = (InitialConditions *)ctx;

//    if (c[0] < initialConditions->length / 2.0) {
//        a_xG[ablate::finiteVolume::CompressibleFlowFields::RHO] = initialConditions->rhoL;
//
//        a_xG[ablate::finiteVolume::CompressibleFlowFields::RHOU + 0] = initialConditions->rhoL * initialConditions->uL;
//
//        PetscReal e = initialConditions->pL / ((initialConditions->gamma - 1.0) * initialConditions->rhoL);
//        PetscReal et = e + 0.5 * PetscSqr(initialConditions->uL);
//        a_xG[ablate::finiteVolume::CompressibleFlowFields::RHOE] = et * initialConditions->rhoL;
//    } else {
//        a_xG[ablate::finiteVolume::CompressibleFlowFields::RHO] = initialConditions->rhoR;
//
//        a_xG[ablate::finiteVolume::CompressibleFlowFields::RHOU + 0] = initialConditions->rhoR * initialConditions->uR;
//
//        PetscReal e = initialConditions->pR / ((initialConditions->gamma - 1.0) * initialConditions->rhoR);
//        PetscReal et = e + 0.5 * PetscSqr(initialConditions->uR);
//        a_xG[ablate::finiteVolume::CompressibleFlowFields::RHOE] = et * initialConditions->rhoR;
//    }
    return 0;
    PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
    // initialize petsc and mpi
    ablate::environment::RunEnvironment::Initialize(&argc, &argv);
    ablate::utilities::PetscUtilities::Initialize();

    {
        // define some initial conditions
        InitialConditions initialConditions{.gamma = 1.4, .Rgas = 287.0, .temperature= 300};

        // setup the run environment
        ablate::parameters::MapParameters runEnvironmentParameters(std::map<std::string, std::string>{{"title", "curvature_AirAir_80x80"}});
        ablate::environment::RunEnvironment::Setup(runEnvironmentParameters);

        auto eos = std::make_shared<ablate::eos::PerfectGas>(std::make_shared<ablate::parameters::MapParameters>(std::map<std::string, std::string>{{"gamma", "1.4"}, {"Rgas", "287.0"}}));

//        auto conservedFieldParameters = std::make_shared<ablate::parameters::Parameters>()

        // determine required fields for finite volume compressible flow
        std::vector<std::shared_ptr<ablate::domain::FieldDescriptor>> fieldDescriptors = {std::make_shared<ablate::finiteVolume::CompressibleFlowFields>(eos)};

//        CompressibleFlowState(std::shared_ptr<ablate::eos::EOS> eos, std::shared_ptr<mathFunctions::MathFunction> temperatureFunction, std::shared_ptr<mathFunctions::MathFunction> pressureFunction,
//                std::shared_ptr<mathFunctions::MathFunction> velocityFunction, std::shared_ptr<mathFunctions::FieldFunction> otherPropertyFunction = {});

//        ablate::finiteVolume::CompressibleFlowFields::CompressibleFlowFields(std::shared_ptr<eos::EOS> eos, std::shared_ptr<domain::Region> region,
//                std::shared_ptr<parameters::Parameters> conservedFieldParameters)

        auto fields = ablate::finiteVolume::fieldFunctions::CompressibleFlowState(std::shared_ptr<ablate::eos::EOS> eos,
                                                                                  std::shared_ptr<ablate::mathFunctions::MathFunction> 300,
                                                                                  std::shared_ptr<ablate::mathFunctions::MathFunction> 10000,
                                                                                  std::shared_ptr<ablate::mathFunctions::MathFunction> 0);

        auto domain =
            std::make_shared<ablate::domain::BoxMesh>("simpleBoxField",
                                                      fieldDescriptors,
                                                      fields,
                                                      std::vector<std::shared_ptr<ablate::domain::modifiers::Modifier>>{std::make_shared<ablate::domain::modifiers::DistributeWithGhostCells>(),
                                                                                                                        std::make_shared<ablate::domain::modifiers::GhostBoundaryCells>()},
                                                      std::vector<int>{40}, //faces
                                                      std::vector<double>{{-2, -2}}, //upper
                                                      std::vector<double>{{2, 2}}, //lower
                                                      std::vector<std::string>{"NONE"} /*boundary*/,
                                                      false /*simplex*/,
                                                      ablate::parameters::MapParameters::Create({{"dm_refine", "1"}}));

        auto serializer =
//        ablate::io::Hdf5MultiFileSerializer::Hdf5MultiFileSerializer(std::shared_ptr<ablate::io::interval::Interval> interval, std::shared_ptr<parameters::Parameters> options)
            std::make_shared<ablate::io::Hdf5MultiFileSerializer>(
                    std::shared_ptr<ablate::io::interval::Interval> 0.1);

        // Set up the flow data
        auto parameters = std::make_shared<ablate::parameters::MapParameters>(std::map<std::string, std::string>{{"cfl", ".4"}});

        // Set the initial conditions for euler
        auto initialCondition = std::make_shared<ablate::mathFunctions::FieldFunction>("euler", ablate::mathFunctions::Create(SetInitialCondition, (void *)&initialConditions));

        // create a time stepper
        auto timeStepper = ablate::solver::TimeStepper(domain,
                                                       serializer,
                                                       ablate::parameters::MapParameters::Create({{"ts_adapt_type", "physicsConstrained"}, {"ts_max_time", "0.5"}, {"ts_dt", "1e-2"}}),
                                                       {},
                                                       std::make_shared<ablate::domain::Initializer>(initialCondition));

        auto boundaryConditions = std::vector<std::shared_ptr<ablate::finiteVolume::boundaryConditions::BoundaryCondition>>{
            std::make_shared<ablate::finiteVolume::boundaryConditions::Ghost>("euler", "wall left", 1, PhysicsBoundary_Euler, (void *)&initialConditions),
            std::make_shared<ablate::finiteVolume::boundaryConditions::Ghost>("euler", "wall right", 2, PhysicsBoundary_Euler, (void *)&initialConditions)};

        // Create a shockTube solver
        auto shockTubeSolver = std::make_shared<ablate::finiteVolume::CompressibleFlowSolver>("compressibleShockTube",
                                                                                              ablate::domain::Region::ENTIREDOMAIN,
                                                                                              nullptr /*options*/,
                                                                                              eos,
                                                                                              parameters,
                                                                                              nullptr /*transportModel*/,
                                                                                              std::make_shared<ablate::finiteVolume::fluxCalculator::Ausm>(),
                                                                                              boundaryConditions /*boundary conditions*/);

        // register the flowSolver with the timeStepper
        timeStepper.Register(
            shockTubeSolver,
            {std::make_shared<ablate::monitors::TimeStepMonitor>(), std::make_shared<ablate::monitors::CurveMonitor>(std::make_shared<ablate::io::interval::FixedInterval>(10), "outputCurve")});

        // Solve the time stepper
        timeStepper.Solve();
    }

    ablate::environment::RunEnvironment::Finalize();
    return 0;
}