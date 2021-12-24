using MultiPrecision;

namespace FresnelIntegral {
    internal static class FindPoint5N8 {
        public static MultiPrecision<Pow2.N8> FresnelC(MultiPrecision<Pow2.N8> x) {
            if (x.Sign == Sign.Minus) {
                return -FresnelC(-x);
            }

            MultiPrecision<Pow2.N8> prev_dx = MultiPrecision<Pow2.N8>.PositiveInfinity;

            while (true) {
                MultiPrecision<Pow2.N8> y = FresnelN8.FresnelC(x) - MultiPrecision<Pow2.N8>.Point5;
                MultiPrecision<Pow2.N8> d = MultiPrecision<Pow2.N8>.CosPI(x * x / 2);
                MultiPrecision<Pow2.N8> dx = y / d;

                if ((dx.Exponent < x.Exponent - MultiPrecision<Pow2.N8>.Bits) ||
                    (MultiPrecision<Pow2.N8>.Abs(prev_dx) <= MultiPrecision<Pow2.N8>.Abs(dx))) {

                    return x;
                }

                x -= dx;
                prev_dx = dx;
            }
        }

        public static MultiPrecision<Pow2.N8> FresnelS(MultiPrecision<Pow2.N8> x) {
            if (x.Sign == Sign.Minus) {
                return -FresnelS(-x);
            }

            MultiPrecision<Pow2.N8> prev_dx = MultiPrecision<Pow2.N8>.PositiveInfinity;

            while (true) {
                MultiPrecision<Pow2.N8> y = FresnelN8.FresnelS(x) - MultiPrecision<Pow2.N8>.Point5;
                MultiPrecision<Pow2.N8> d = MultiPrecision<Pow2.N8>.SinPI(x * x / 2);
                MultiPrecision<Pow2.N8> dx = y / d;

                if ((dx.Exponent < x.Exponent - MultiPrecision<Pow2.N8>.Bits) ||
                    (MultiPrecision<Pow2.N8>.Abs(prev_dx) <= MultiPrecision<Pow2.N8>.Abs(dx))) {

                    return x;
                }

                x -= dx;
                prev_dx = dx;
            }
        }
    }
}
