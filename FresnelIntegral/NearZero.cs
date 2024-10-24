using MultiPrecision;

namespace FresnelIntegral {
    internal static class NearZero<N, M> where N : struct, IConstant where M : struct, IConstant {
        public static MultiPrecision<N> FresnelC(MultiPrecision<N> x, int max_terms = 4096) {
            if (x.Sign == Sign.Minus) {
                return -FresnelC(-x);
            }
            if (MultiPrecision<N>.IsZero(x)) {
                return 0;
            }

            MultiPrecision<M> x_ex = x.Convert<M>();

            MultiPrecision<M> v = x_ex * x_ex * MultiPrecision<M>.PI;
            MultiPrecision<M> v2 = v * v, v4 = v2 * v2;

            MultiPrecision<M> s = 0, u = x_ex;

            long k = 0;

            for (long conv_times = 0; k < max_terms && conv_times < 3; k++) {
                MultiPrecision<M> f = MultiPrecision<M>.Ldexp(u * TaylorSeries<M>.Value(checked((int)(4 * k))), -4 * k);
                MultiPrecision<M> ds = f * (MultiPrecision<M>.Div(1, 8 * k + 1) - v2 / checked(4 * (8 * k + 5) * (4 * k + 1) * (4 * k + 2)));

                if ((x <= 1 && ds.Exponent < s.Exponent - MultiPrecision<N>.Bits) ||
                    (x > 1 && ds.Exponent < -MultiPrecision<N>.Bits - 1)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                s += ds;
                u *= v4;

                if (s.Exponent > MultiPrecision<M>.Bits - MultiPrecision<N>.Bits) {
                    return MultiPrecision<N>.NaN;
                }
            }

            if (k < max_terms) {
                return s.Convert<N>();
            }
            else {
                return MultiPrecision<N>.NaN;
            }
        }

        public static MultiPrecision<N> FresnelS(MultiPrecision<N> x, int max_terms = 1024) {
            if (x.Sign == Sign.Minus) {
                return -FresnelS(-x);
            }
            if (MultiPrecision<N>.IsZero(x)) {
                return 0;
            }

            MultiPrecision<M> x_ex = x.Convert<M>();

            MultiPrecision<M> v = x_ex * x_ex * MultiPrecision<M>.PI;
            MultiPrecision<M> v2 = v * v, v4 = v2 * v2;

            MultiPrecision<M> s = 0, u = v * x_ex / 2;

            long k = 0;

            for (long conv_times = 0; k < max_terms && conv_times < 3; k++) {
                MultiPrecision<M> f = MultiPrecision<M>.Ldexp(u * TaylorSeries<M>.Value(checked((int)(4 * k + 1))), -4 * k);
                MultiPrecision<M> ds = f * (MultiPrecision<M>.Div(1, 8 * k + 3) - v2 / checked(4 * (8 * k + 7) * (4 * k + 2) * (4 * k + 3)));

                if ((x <= 1 && ds.Exponent < s.Exponent - MultiPrecision<N>.Bits) ||
                    (x > 1 && ds.Exponent < -MultiPrecision<N>.Bits - 1)) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                s += ds;
                u *= v4;

                if (s.Exponent > MultiPrecision<M>.Bits - MultiPrecision<N>.Bits) {
                    return MultiPrecision<N>.NaN;
                }
            }

            if (k < max_terms) {
                return s.Convert<N>();
            }
            else {
                return MultiPrecision<N>.NaN;
            }
        }
    }
}
