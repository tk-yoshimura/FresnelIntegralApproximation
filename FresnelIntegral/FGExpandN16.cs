using System;
using MultiPrecision;

namespace FresnelIntegral {
    internal static class FGExpandN16 {
        private const double threshold = 15.75d;

        public static (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) Value(MultiPrecision<Pow2.N16> x) {
            if (x.Sign == Sign.Minus) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            if (x < threshold) {
                MultiPrecision<Pow2.N16> fc = NearZero<Pow2.N16, Pow2.N64>.FresnelC(x);
                MultiPrecision<Pow2.N16> fs = NearZero<Pow2.N16, Pow2.N64>.FresnelS(x);

                MultiPrecision<Pow2.N16> theta = x * x / 2;
                MultiPrecision<Pow2.N16> cos = MultiPrecision<Pow2.N16>.CosPI(theta);
                MultiPrecision<Pow2.N16> sin = MultiPrecision<Pow2.N16>.SinPI(theta);

                MultiPrecision<Pow2.N16> f = ((1 - 2 * fs) * cos - (1 - 2 * fc) * sin) / 2;
                MultiPrecision<Pow2.N16> g = ((1 - 2 * fs) * sin + (1 - 2 * fc) * cos) / 2;

                return (f, g);
            }
            else {
                return Limit<Pow2.N16>.Coef(x);
            }
        }
    }
}
