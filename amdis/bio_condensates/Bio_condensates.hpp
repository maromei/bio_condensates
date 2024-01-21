#pragma once

namespace AMDiS::Bio_condensates {

    template <typename M>
    auto doubleWell(M phi) {
        return 18 * pow<2>(phi) * pow<2>(1. - phi);
    }

    template <typename M>
    auto doubleWellDeriv(M phi) {
        return 36. * (2*pow<3>(phi) - 3*pow<2>(phi) + phi);
    }

    template <typename M>
    auto doubleWellSecondDeriv(M phi) {
        return 36. * (6*pow<2>(phi) - 6*phi + 1.);
    }

};
