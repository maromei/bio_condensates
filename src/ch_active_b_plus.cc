#include <string>
#include <iostream>

#include <amdis/AMDiS.hpp>
#include <amdis/AdaptInstationary.hpp>
#include <amdis/LocalOperators.hpp>
#include <amdis/ProblemInstat.hpp>
#include <amdis/ProblemStat.hpp>
#include <amdis/GridFunctions.hpp>
#include <amdis/Marker.hpp>

#include <dune/grid/albertagrid.hh>

#include <amdis/bio_condensates/Bio_condensates.hpp>

using namespace AMDiS;
using namespace Bio_condensates;

using Grid = Dune::AlbertaGrid<GRIDDIM, WORLDDIM>;
using Param = LagrangeBasis<Grid, 1, 1, 1, 1, 1, 1>;

int main(int argc, char** argv) {

    Environment env(argc, argv);

    ProblemStat<Param> prob("ch");
    prob.initialize(INIT_ALL);

    ProblemInstat<Param> probInstat("ch", prob);
    probInstat.initialize(INIT_UH_OLD);

    AdaptInfo adaptInfo("adapt");

    /***************
    ** Parameters **
    ***************/

    auto invTau = std::ref(probInstat.invTau());

    auto phi = prob.solution(0);
    auto phiOld = probInstat.oldSolution(0);

    auto K = prob.solution(1);
    auto mu = prob.solution(2);
    auto D = prob.solution(3);
    auto Dx = prob.solution(4);
    auto Dy = prob.solution(5);

    double mobility = Parameters::get<double>("parameters->mobility").value_or(1.0);
    double eps = Parameters::get<double>("parameters->epsilon").value_or(0.02);
    double lambda = Parameters::get<double>("parameters->lambda").value_or(1.0);
    double gamma = Parameters::get<double>("parameters->gamma").value_or(1.0);
    double avg_density = Parameters::get<double>("parameters->avg_density").value_or(1.0);
    std::string initType = Parameters::get<std::string>("parameters->initType").value_or("singleBubble");

    /**************
    ** Operators **
    **************/

    // ##############
    // ### Matrix ###
    // ##############
    // ### M_00 --> eq1, phi
    auto opDnx = makeOperator(
        tag::test_partialtrial{ 0 },
        - gamma * Dx
    );
    auto opDny = makeOperator(
        tag::test_partialtrial{ 1 },
        - gamma * Dy
    );

    prob.addMatrixOperator(zot(invTau), 0, 0);
    prob.addMatrixOperator(opDnx, 0, 0);
    prob.addMatrixOperator(opDny, 0, 0);

    // ### M_01 --> eq1, K
    prob.addMatrixOperator(sot(1.0), 0, 1);

    // ## M_0 Zeros
    prob.addMatrixOperator(zot(0.0), 0, 2);
    prob.addMatrixOperator(zot(0.0), 0, 3);
    prob.addMatrixOperator(zot(0.0), 0, 4);
    prob.addMatrixOperator(zot(0.0), 0, 5);

    // ###########################################

    // ### M_10 --> eq2, phi
    prob.addMatrixOperator(zot(
        - lambda / eps * doubleWellDeriv(phi)
    ), 1, 0);

    // ### M_11 --> eq2, K
    prob.addMatrixOperator(zot(1.0), 1, 1);

    // ### M_12 --> eq2, mu
    prob.addMatrixOperator(zot(-mobility), 1, 2);

    // ### M_1 Zeros
    prob.addMatrixOperator(zot(0.0), 1, 3);
    prob.addMatrixOperator(zot(0.0), 1, 4);
    prob.addMatrixOperator(zot(0.0), 1, 5);

    // ###########################################

    // ### M_20 --> eq3, phi
    prob.addMatrixOperator(zot(
        - 1. / eps * doubleWellSecondDeriv(phi)
    ), 2, 0);

    // ### M_22 --> eq3, mu
    prob.addMatrixOperator(zot(1.0), 2, 2);

    // ### M_23 --> eq3, D
    prob.addMatrixOperator(zot(eps), 2, 3);

    // ### M_2 Zeros
    prob.addMatrixOperator(zot(0.0), 2, 1);
    prob.addMatrixOperator(zot(0.0), 2, 4);
    prob.addMatrixOperator(zot(0.0), 2, 5);

    // ###########################################

    // ### M_30 --> eq4, phi
    prob.addMatrixOperator(sot(1.0), 3, 0);

    // ### M_33 --> eq4, D
    prob.addMatrixOperator(zot(1.0), 3, 3);

    // ### M_3 Zeros
    prob.addMatrixOperator(zot(0.0), 3, 1);
    prob.addMatrixOperator(zot(0.0), 3, 2);
    prob.addMatrixOperator(zot(0.0), 3, 4);
    prob.addMatrixOperator(zot(0.0), 3, 5);

    // ###########################################

    // ### M_43 --> eq5, D
    auto opDpartialx = makeOperator(
        tag::test_partialtrial{ 0 },
        1.0
    );
    prob.addMatrixOperator(opDpartialx, 4, 3);

    // ### M_44 --> eq5, Dx
    prob.addMatrixOperator(zot(1.0), 4, 4);

    // ### M_4 Zeros
    prob.addMatrixOperator(zot(0.0), 4, 0);
    prob.addMatrixOperator(zot(0.0), 4, 1);
    prob.addMatrixOperator(zot(0.0), 4, 2);
    prob.addMatrixOperator(zot(0.0), 4, 5);

    // ###########################################

    // ### M_53 --> eq6, D
    auto opDpartialy = makeOperator(
        tag::test_partialtrial{ 1 },
        1.0
    );
    prob.addMatrixOperator(opDpartialy, 5, 3);

    // ### M_55 --> eq6, Dy
    prob.addMatrixOperator(zot(1.0), 5, 5);

    // ### M_5 Zeros
    prob.addMatrixOperator(zot(0.0), 5, 0);
    prob.addMatrixOperator(zot(0.0), 5, 1);
    prob.addMatrixOperator(zot(0.0), 5, 2);
    prob.addMatrixOperator(zot(0.0), 5, 4);

    // ##############
    // ### Vector ###
    // ##############
    // ### f_0
    prob.addVectorOperator(zot(phiOld * invTau), 0);
    prob.addVectorOperator(zot(- gamma * pow<2>(D)), 0);

    // ### f_1
    prob.addVectorOperator(zot(
        lambda * (doubleWell(phi) - doubleWellDeriv(phi) * phi)
    ), 1);

    // ### f_2
    prob.addVectorOperator(zot(
        1./eps * (doubleWellDeriv(phi) - doubleWellSecondDeriv(phi) * phi)
    ), 2);

    // ### f Zeros
    prob.addVectorOperator(zot(0.0), 3);
    prob.addVectorOperator(zot(0.0), 4);
    prob.addVectorOperator(zot(0.0), 5);


    /************************
    ** Boundary Management **
    ************************/

    /*
    // dirichlet boundary condition
    auto predicate = [](auto const& x){
        return x[0] < 1.e-8 || x[1] < 1.e-8;
    };
    auto dbcValues = [](auto const& x){ return 0.0; }; // set value
    prob.addDirichletBC(predicate, 1, 1, dbcValues);
    */

    /*******************
    ** Grid refinment **
    *******************/

    // Define an indicator for local adaptive refinement
    int ref_int  = Parameters::get<int>("refinement->interface").value_or(10);
    int ref_bulk = Parameters::get<int>("refinement->bulk").value_or(2);
    int ref_init = Parameters::get<int>("refinement->init").value_or(
        ref_int - ref_bulk
    );
    auto indicator = [ref_int, ref_bulk](double phi) {
        return phi > 0.05 && phi < 0.95 ? ref_int : ref_bulk;
    };

    // Construct a grid marker for the automatic refinement
    GridFunctionMarker marker("interface", prob.grid(), invokeAtQP(indicator,phi));
    prob.addMarker(marker);

    /**********************
    ** Initial Condition **
    **********************/

    // Define an initial condition for the phase field

    if (initType == "uniformAvgDensity")
        eps = avg_density;

    // Perform an initial refinement
    for (int i = 0; i < ref_init; ++i) {
        Initializer::initilizeWithType(initType, eps, phi);
        prob.markElements(adaptInfo);
        prob.adaptGrid(adaptInfo);
    }

    /*******************
    ** Run Simulation **
    *******************/

    AdaptInstationary adapt("adapt", prob, adaptInfo, probInstat, adaptInfo);
    adapt.adapt();

    return 0;
}
