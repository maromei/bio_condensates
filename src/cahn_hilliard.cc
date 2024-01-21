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

    const int bubbleNum = 100;
    double radiiCoords[bubbleNum][3] = {
        {0.5769, 1.8439, -5.8772},
        {0.1409, 2.6895, 6.4982},
        {0.6307, -2.1079, -2.2185},
        {0.5924, 0.9430, 3.0513},
        {0.7654, 0.4854, -8.1086},
        {0.1106, -5.1587, 4.9410},
        {0.1319, 2.3709, 9.3420},
        {0.7661, 2.9311, -9.1401},
        {0.1489, -0.2121, 1.3676},
        {0.5212, 5.4661, -5.1190},
        {0.4646, -6.6962, 7.8503},
        {0.1150, -1.3441, 2.1413},
        {0.6670, -5.6697, 0.0136},
        {0.8633, 2.4743, 2.9555},
        {0.2365, 8.1286, 3.3518},
        {0.5132, -8.4171, 3.5795},
        {0.3782, -1.7941, 9.0278},
        {0.8907, 6.4068, 0.2169},
        {0.1946, -6.8398, -7.7842},
        {0.7346, -2.1301, 0.9893},
        {0.5153, -4.4434, 0.9707},
        {0.2546, -8.7075, -1.5459},
        {0.3194, -7.7492, 0.7388},
        {0.2894, -6.6017, 6.0109},
        {0.5235, -2.8327, -8.0448},
        {0.7373, -7.7991, 6.5588},
        {0.5810, -7.1702, 3.2402},
        {0.3731, -9.0634, -5.8855},
        {0.1201, -8.4406, -2.1097},
        {0.5502, -9.1736, 1.3434},
        {0.9296, -3.2818, 7.2582},
        {0.5708, 6.7219, 9.2151},
        {0.6584, 5.4203, 7.1790},
        {0.9913, 3.6965, -0.9198},
        {0.7167, 7.3709, -2.7517},
        {0.2323, 4.3795, 0.4904},
        {0.4985, 8.7023, 1.7271},
        {0.6033, -0.8729, -4.8068},
        {0.2810, 8.0136, 8.4649},
        {0.3518, -5.7571, -7.8699},
        {0.8752, 8.5418, 5.0438},
        {0.2189, 9.1332, -1.4379},
        {0.1044, 6.1381, -6.8576},
        {0.3374, -0.6175, 0.1220},
        {0.1000, -4.7168, 3.2567},
        {0.5374, -1.1763, -8.4539},
        {0.3305, -4.7043, -2.7440},
        {0.5108, 6.3462, 5.9681},
        {0.5937, 6.8513, -8.5022},
        {0.3253, -6.9187, -5.6386},
        {0.9534, -1.4857, 3.4801},
        {0.3015, -8.3452, -0.9027},
        {0.2892, -4.6643, -5.1180},
        {0.7940, 5.3273, -3.3738},
        {0.5819, 8.4542, 0.5429},
        {0.5969, 8.6153, -5.1485},
        {0.1216, -3.9491, 2.3190},
        {0.4489, -5.8995, -9.0152},
        {0.4884, 9.0267, -6.6658},
        {0.4493, -5.4026, -5.9511},
        {0.1110, 5.0842, -6.9348},
        {0.1299, 4.4878, 0.9859},
        {0.7216, -5.9919, -1.7508},
        {0.3497, -1.7721, 7.4068},
        {0.4971, 3.0914, 9.4610},
        {0.8423, 6.3777, 3.6327},
        {0.2012, 0.6603, -0.8965},
        {0.7319, 2.8380, 7.9073},
        {0.1919, -0.4848, 7.8054},
        {0.1321, 2.1988, -0.9064},
        {0.2324, 4.5766, -8.9738},
        {0.3027, -0.7355, -1.2029},
        {0.5282, -9.3190, -3.6901},
        {0.2309, 5.0897, 9.1607},
        {0.8542, 0.4883, 6.6839},
        {0.3201, -3.2233, -2.2205},
        {0.1794, 3.6977, -4.1308},
        {0.5936, 2.9422, 5.0986},
        {0.2819, 2.3848, -6.6163},
        {0.5569, 7.7972, -0.9753},
        {0.1306, -6.2841, 9.6538},
        {0.4295, -1.6425, -7.1125},
        {0.2171, 1.4933, 4.2553},
        {0.1200, 7.3516, 1.1898},
        {0.3843, -9.0786, -7.9037},
        {0.2027, -8.2572, -6.7324},
        {0.2702, -9.4221, 2.7104},
        {0.1167, -7.5730, -6.6561},
        {0.3000, -8.0532, 4.9574},
        {0.2975, -2.1508, -4.9519},
        {0.4531, -9.3985, -6.7775},
        {0.5783, -3.8256, 4.8456},
        {0.5635, -6.4369, 1.0245},
        {0.8317, -3.9839, -6.2090},
        {0.6426, -5.2679, 6.6052},
        {0.1979, 9.6136, -2.3107},
        {0.2472, 4.7564, 3.6202},
        {0.3165, 9.2320, -9.5297},
        {0.1231, -3.5481, -7.9925},
        {0.3051, -7.5347, -8.1407}
    };

    auto phase = [eps](
        auto const& worldVec,
        double radius,
        double offset_x,
        double offset_y
    ) {

        auto x = worldVec[0] - offset_x;
        auto y = worldVec[1] - offset_y;
        auto dist = sqrt(x*x + y*y) - radius;

        return 0.5*(1 + std::tanh(-3.0 * dist/eps));
    };

    auto initPhi = [phase, radiiCoords, bubbleNum](auto const& x) {
        auto phi = phase(x, radiiCoords[0][0], radiiCoords[0][1], radiiCoords[0][2]);
        for (int i = 1; i < bubbleNum; i++) {
            phi += phase(x, radiiCoords[i][0], radiiCoords[i][1], radiiCoords[i][2]);
        }
        return phi;
    };

    auto defaultInitCircle = [phase](auto const& x) {
        return phase(x, 0.5, 0., 0.);
    };

    // Perform an initial refinement
    for (int i = 0; i < ref_init; ++i) {
        phi << initPhi;
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
