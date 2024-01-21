#pragma once

#include <string>
#include <stdexcept>
#include <functional>

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

    namespace Initializer {

        const int bubbleNum = 100;
        double radiiCoordsBigBubble[bubbleNum][3] = {
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

        template <typename T>
        auto phase(
            T const& worldVec,
            double radius,
            double offset_x,
            double offset_y,
            double eps
        ) {

            auto x = worldVec[0] - offset_x;
            auto y = worldVec[1] - offset_y;
            auto dist = sqrt(x*x + y*y) - radius;

            return 0.5*(1 + std::tanh(-3.0 * dist/eps));
        };

        auto defaultInitCircle(double eps) {

            auto initCircleLambda = [eps](auto const& x) {
                return Initializer::phase(x, 5, 0., 0., eps);
            };

            return initCircleLambda;
        }

        auto bigBubble(double eps) {

            auto initPhi = [eps](auto const& x) {
                auto phi = Initializer::phase(
                    x,
                    Initializer::radiiCoordsBigBubble[0][0],
                    Initializer::radiiCoordsBigBubble[0][1],
                    Initializer::radiiCoordsBigBubble[0][2],
                    eps
                );
                for (int i = 1; i < Initializer::bubbleNum; i++) {
                    phi += phase(
                        x,
                        Initializer::radiiCoordsBigBubble[i][0],
                        Initializer::radiiCoordsBigBubble[i][1],
                        Initializer::radiiCoordsBigBubble[i][2],
                        eps
                    );
                }
                return phi;
            };

            return initPhi;

        }

        auto uniformAvgDensity(double avgDensity) {

            auto init = [avgDensity](auto const& x) {
                return avgDensity;
            };

            return init;
        }

        template <typename T>
        void initilizeWithType(std::string type, double eps, T phi) {

            if (type == "manyBubbles") {
                phi << Initializer::bigBubble(eps);
            } else if (type == "singleBubble") {
                phi << Initializer::defaultInitCircle(eps);
            } else if (type == "uniformAvgDensity") {
                phi << Initializer::uniformAvgDensity(eps);
            } else {
                throw std::runtime_error("Invalid init type.");
            }

            return;
        }

    };

};
