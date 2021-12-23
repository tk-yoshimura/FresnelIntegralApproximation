using MultiPrecision;

namespace FresnelIntegral {
    internal static class Limit<N> where N : struct, IConstant {

        public static (MultiPrecision<N> c, MultiPrecision<N> s) Fresnel(MultiPrecision<N> x, int max_terms = 64) {
            (MultiPrecision<N> p, MultiPrecision<N> q) = Coef(x, max_terms);

            MultiPrecision<N> theta = x * x / 2;
            MultiPrecision<N> cos = MultiPrecision<N>.CosPI(theta);
            MultiPrecision<N> sin = MultiPrecision<N>.SinPI(theta);

            MultiPrecision<N> c = MultiPrecision<N>.Point5 + sin * p - cos * q;
            MultiPrecision<N> s = MultiPrecision<N>.Point5 - cos * p - sin * q;

            return (c, s);
        }

        public static (MultiPrecision<N> p, MultiPrecision<N> q) Coef(MultiPrecision<N> x, int max_terms = 64) {
            MultiPrecision<N> x_ex = x.Convert<N>();

            MultiPrecision<N> v = MultiPrecision<N>.Rcp(x_ex * x_ex * MultiPrecision<N>.PI);
            MultiPrecision<N> v2 = v * v, v4 = v2 * v2;

            MultiPrecision<N> p = 0, q = 0;
            MultiPrecision<N> a = MultiPrecision<N>.Rcp(x_ex * MultiPrecision<N>.PI);
            MultiPrecision<N> b = MultiPrecision<N>.Rcp(x_ex * x_ex * x_ex * MultiPrecision<N>.PI * MultiPrecision<N>.PI);

            for (long k = 0; k < max_terms; k++) {
                MultiPrecision<N> s = checked((8 * k + 1) * (8 * k + 3)) * v2;
                MultiPrecision<N> t = checked((8 * k + 3) * (8 * k + 5)) * v2;

                if (s > 1 || t > 1) {
                    break;
                }

                MultiPrecision<N> dp = a * (1 - s) * RSeries<N>.Value(checked((int)(4 * k)));
                MultiPrecision<N> dq = b * (1 - t) * RSeries<N>.Value(checked((int)(4 * k + 1)));

                p += dp;
                q += dq;
                a *= v4;
                b *= v4;

                if ((dp.Exponent < p.Exponent - MultiPrecision<N>.Bits) &&
                    (dq.Exponent < q.Exponent - MultiPrecision<N>.Bits)) {
                    return (p, q);
                }
            }

            return (MultiPrecision<N>.NaN, MultiPrecision<N>.NaN);
        }
    }
}
