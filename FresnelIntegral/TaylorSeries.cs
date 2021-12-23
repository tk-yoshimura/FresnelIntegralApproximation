using MultiPrecision;
using System;
using System.Collections.Generic;

namespace FresnelIntegral {
    internal static class TaylorSeries<N> where N: struct, IConstant {
        private static MultiPrecision<N> v;
        private readonly static List<MultiPrecision<N>> table;

        static TaylorSeries() {
            v = 1;
            table = new List<MultiPrecision<N>>(new MultiPrecision<N>[] { 1, 1 });
        }

        public static MultiPrecision<N> Value(int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            if (n < table.Count) {
                return table[n];
            }

            for (int m = table.Count; m <= n; m++) {
                v *= m;                
                table.Add(MultiPrecision<N>.Rcp(v));

                if (table[^1].IsZero) {
                    throw new ArithmeticException("taylor series is underflow.");
                }
            }

            return table[n];
        }
    }
}
