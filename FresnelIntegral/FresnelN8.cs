using MultiPrecision;

namespace FresnelIntegral {
    internal static class FresnelN8 {
        private const double threshold = 10.75d;

        public static (MultiPrecision<Pow2.N8> cos, MultiPrecision<Pow2.N8> sin) Value(MultiPrecision<Pow2.N8> x) {
            if (x.Sign == Sign.Minus) {
                (MultiPrecision<Pow2.N8> cos, MultiPrecision<Pow2.N8> sin) = Value(-x);

                return (-cos, -sin);
            }

            if (x < threshold) {
                MultiPrecision<Pow2.N8> cos = NearZero<Pow2.N8, Pow2.N16>.FresnelC(x);
                MultiPrecision<Pow2.N8> sin = NearZero<Pow2.N8, Pow2.N16>.FresnelS(x);

                return (cos, sin);
            }
            else {
                return Limit<Pow2.N8>.Fresnel(x);
            }
        }

        public static MultiPrecision<Pow2.N8> FresnelC(MultiPrecision<Pow2.N8> x) {
            if (x.Sign == Sign.Minus) {
                return -FresnelC(-x);
            }

            if (x < threshold) {
                return NearZero<Pow2.N8, Pow2.N16>.FresnelC(x);
            }
            else {
                return Limit<Pow2.N8>.Fresnel(x).c;
            }
        }

        public static MultiPrecision<Pow2.N8> FresnelS(MultiPrecision<Pow2.N8> x) {
            if (x.Sign == Sign.Minus) {
                return -FresnelS(-x);
            }

            if (x < threshold) {
                return NearZero<Pow2.N8, Pow2.N16>.FresnelS(x);
            }
            else {
                return Limit<Pow2.N8>.Fresnel(x).s;
            }
        }
    }
}
