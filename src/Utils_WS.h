#ifndef UTILS_WS_H
#define UTILS_WS_H

namespace Utils_WS
{
    constexpr float lumiRelError = 0.017;

    /**
     * @brief Calculate significance given p-value
     *
     *             -1
     * Z = √2 x erf  (1 - 2p)
     *
     */
    double pvalueToSignificance(double p);

    /**
     * @brief Calculate p-value given significance
     *
     *     1               Z
     * p = - x [ 1 - erf (---)]
     *     2              √2
     *
     */
    double significanceToPvalue(double z);

    /**
     * @todo Check the definition
     *
     */
    double q0ToSignificance(double q0);

    /**
     * @todo Check the definition
     *
     */
    double nllToQ0(double nllMu0, double nllMuHat, double muHat);

}

#endif
