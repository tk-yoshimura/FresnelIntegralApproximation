using MultiPrecision;
using System;
using System.Collections.Generic;

namespace FresnelIntegral {
    internal static class FGDiffSeries<N> where N : struct, IConstant {
        public static (MultiPrecision<N>[] fs, MultiPrecision<N>[] gs) Diff(MultiPrecision<N> x, MultiPrecision<N> f, MultiPrecision<N> g, int n) {
            ArgumentOutOfRangeException.ThrowIfLessThan(n, 2);

            List<MultiPrecision<N>> fs = [f, -MultiPrecision<N>.PI * x * g];
            List<MultiPrecision<N>> gs = [g, +MultiPrecision<N>.PI * x * f - 1];

            for (int i = 2; i <= n; i++) {
                MultiPrecision<N> df = -MultiPrecision<N>.PI * ((i - 1) * gs[i - 2] + x * gs[i - 1]);
                MultiPrecision<N> dg = +MultiPrecision<N>.PI * ((i - 1) * fs[i - 2] + x * fs[i - 1]);

                fs.Add(df);
                gs.Add(dg);
            }

            return (fs.ToArray(), gs.ToArray());
        }
    }
}
