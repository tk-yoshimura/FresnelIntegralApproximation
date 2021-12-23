using MultiPrecision;

namespace FresnelIntegral {
    internal static class FresnelC<N, M> where N : struct, IConstant where M: struct, IConstant {
        public static MultiPrecision<N> Value(MultiPrecision<N> x, int max_terms = 1024) {
            if (x.Sign == Sign.Minus) {
                return -Value(-x);
            }
            if (x.IsZero) {
                return 0;
            }

            MultiPrecision<M> x_ex = x.Convert<M>();

            MultiPrecision<M> v = x_ex * x_ex * MultiPrecision<M>.PI;
            MultiPrecision<M> v2 = v * v, v4 = v2 * v2;

            MultiPrecision<M> s = 0, u = x_ex;

            int k = 0;
            
            for (int conv_times = 0; k < max_terms && conv_times < 3; k++) {
                MultiPrecision<M> f = MultiPrecision<M>.Ldexp(u * TaylorSeries<M>.Value(4 * k), -4 * k);
                MultiPrecision<M> ds = f * (MultiPrecision<M>.Div(1, 8 * k + 1) - v2 / (4 * (8 * k + 5) * (4 * k + 1) * (4 * k + 2)));

                if (ds.Exponent < -MultiPrecision<N>.Bits - 1) {
                    conv_times++;
                }
                else {
                    conv_times = 0;
                }

                s += ds;
                u *= v4;

                if (s.Exponent > MultiPrecision<N>.Bits) { 
                    return MultiPrecision<N>.NaN;
                }
            }

            var n = TaylorSeries<M>.Value(4 * max_terms);

            if (k < max_terms) {
                return s.Convert<N>();
            }
            else {
                return MultiPrecision<N>.NaN;
            }
        }
    }
}
