using MultiPrecision;
using System;
using System.Collections.Generic;
using System.Numerics;

namespace FresnelIntegral {
    internal static class TaylorSeries<N> where N : struct, IConstant {
        private static BigInteger v;
        private static readonly List<MultiPrecision<N>> table;

        static TaylorSeries() {
            v = 1;
            table = new List<MultiPrecision<N>>([1, 1]);
        }

        public static MultiPrecision<N> Value(int n) {
            ArgumentOutOfRangeException.ThrowIfNegative(n);

            if (n < table.Count) {
                return table[n];
            }

            for (int m = table.Count; m <= n; m++) {
                v *= m;
                table.Add(MultiPrecision<N>.Rcp(v));

                if (MultiPrecision<N>.IsZero(table[^1])) {
                    throw new ArithmeticException("taylor series is underflow.");
                }
            }

            return table[n];
        }
    }
}
