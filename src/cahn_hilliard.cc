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
using Param = LagrangeBasis<Grid, 1, 1>;

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

    double mobility = Parameters::get<double>("parameters->mobility").value_or(1.0);
    double eps = Parameters::get<double>("parameters->epsilon").value_or(0.02);
    std::string initType = Parameters::get<std::string>("parameters->initType").value_or("singleBubble");

    /**************
    ** Operators **
    **************/

    // ### Matrix ###
    // ### M_00

    prob.addMatrixOperator(zot(invTau), 0, 0);

    // ### M_01

    prob.addMatrixOperator(sot(mobility), 0, 1);

    // ### M_10

    prob.addMatrixOperator(sot(-eps), 1, 0);
    auto linDoubleWellOperator = zot(
        - 1./eps * doubleWellSecondDeriv(phi)
    );
    prob.addMatrixOperator(linDoubleWellOperator, 1, 0);

    // ### M_11

    prob.addMatrixOperator(zot(1.0), 1, 1);

    // ### Vector ###
    // ### f_0

    prob.addVectorOperator(zot(phiOld * invTau), 0);

    // ### f_1

    auto nonLinDoubleWellOperator = zot(
        1./eps * doubleWellDeriv(phi) -
        1./eps * doubleWellSecondDeriv(phi) * phi
    );
    prob.addVectorOperator(nonLinDoubleWellOperator, 1);

    /************************
    ** Boundary Management **
    ************************/

    //prob.boundaryManager()->setBoxBoundary()

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
