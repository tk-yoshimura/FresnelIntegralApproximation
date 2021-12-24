using System;
using MultiPrecision;

namespace FresnelIntegral {
    internal static class FGExpandN8 {
        private const double threshold = 10.75d;

        public static (MultiPrecision<Pow2.N8> f, MultiPrecision<Pow2.N8> g) Value(MultiPrecision<Pow2.N8> x) {
            if (x.Sign == Sign.Minus) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            if (x < threshold) {
                MultiPrecision<Pow2.N8> fc = NearZero<Pow2.N8, Pow2.N16>.FresnelC(x);
                MultiPrecision<Pow2.N8> fs = NearZero<Pow2.N8, Pow2.N16>.FresnelS(x);

                MultiPrecision<Pow2.N8> theta = x * x / 2;
                MultiPrecision<Pow2.N8> cos = MultiPrecision<Pow2.N8>.CosPI(theta);
                MultiPrecision<Pow2.N8> sin = MultiPrecision<Pow2.N8>.SinPI(theta);

                MultiPrecision<Pow2.N8> f = ((1 - 2 * fs) * cos - (1 - 2 * fc) * sin) / 2;
                MultiPrecision<Pow2.N8> g = ((1 - 2 * fs) * sin + (1 - 2 * fc) * cos) / 2;

                return (f, g);
            }
            else {
                return Limit<Pow2.N8>.Coef(x);
            }
        }
    }
}
