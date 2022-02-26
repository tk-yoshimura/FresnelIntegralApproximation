using MultiPrecision;
using System;

namespace FresnelIntegral {
    internal static class FGExpandN32 {
        private const double threshold = 32;

        public static (MultiPrecision<Pow2.N32> f, MultiPrecision<Pow2.N32> g) Value(MultiPrecision<Pow2.N32> x) {
            if (x.Sign == Sign.Minus) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            if (x <= threshold) {
                MultiPrecision<Pow2.N32> fc = NearZero<Pow2.N32, Pow2.N128>.FresnelC(x, max_terms: 2048);
                MultiPrecision<Pow2.N32> fs = NearZero<Pow2.N32, Pow2.N128>.FresnelS(x, max_terms: 2048);

                MultiPrecision<Pow2.N32> theta = x * x / 2;
                MultiPrecision<Pow2.N32> cos = MultiPrecision<Pow2.N32>.CosPI(theta);
                MultiPrecision<Pow2.N32> sin = MultiPrecision<Pow2.N32>.SinPI(theta);

                MultiPrecision<Pow2.N32> f = ((1 - 2 * fs) * cos - (1 - 2 * fc) * sin) / 2;
                MultiPrecision<Pow2.N32> g = ((1 - 2 * fs) * sin + (1 - 2 * fc) * cos) / 2;

                return (f, g);
            }
            else {
                return Limit<Pow2.N32>.Coef(x);
            }
        }
    }
}
