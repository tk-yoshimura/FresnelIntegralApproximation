using MultiPrecision;

namespace FresnelIntegral {
    internal static class FresnelN16 {
        private const double threshold = 15.75d;

        public static (MultiPrecision<Pow2.N16> cos, MultiPrecision<Pow2.N16> sin) Value(MultiPrecision<Pow2.N16> x) {
            if (x.Sign == Sign.Minus) {
                (MultiPrecision<Pow2.N16> cos, MultiPrecision<Pow2.N16> sin) = Value(-x);

                return (-cos, -sin);
            }

            if (x < threshold) {
                MultiPrecision<Pow2.N16> cos = NearZero<Pow2.N16, Pow2.N64>.FresnelC(x);
                MultiPrecision<Pow2.N16> sin = NearZero<Pow2.N16, Pow2.N64>.FresnelS(x);

                return (cos, sin);
            }
            else {
                return Limit<Pow2.N16>.Fresnel(x);
            }
        }

        public static MultiPrecision<Pow2.N16> FresnelC(MultiPrecision<Pow2.N16> x) {
            if (x.Sign == Sign.Minus) {
                return -FresnelC(-x);
            }

            if (x < threshold) {
                return NearZero<Pow2.N16, Pow2.N64>.FresnelC(x);
            }
            else {
                return Limit<Pow2.N16>.Fresnel(x).c;
            }
        }

        public static MultiPrecision<Pow2.N16> FresnelS(MultiPrecision<Pow2.N16> x) {
            if (x.Sign == Sign.Minus) {
                return -FresnelS(-x);
            }

            if (x < threshold) {
                return NearZero<Pow2.N16, Pow2.N64>.FresnelS(x);
            }
            else {
                return Limit<Pow2.N16>.Fresnel(x).s;
            }
        }
    }
}
