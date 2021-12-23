using MultiPrecision;
using System;
using System.Collections.Generic;
using System.Numerics;

namespace FresnelIntegral {
    internal class RSeries<N> where N : struct, IConstant {
        private static BigInteger v;
        private static readonly List<MultiPrecision<N>> table;

        static RSeries() {
            v = 1;
            table = new List<MultiPrecision<N>>(new MultiPrecision<N>[] { 1 });
        }

        public static MultiPrecision<N> Value(int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            if (n < table.Count) {
                return table[n];
            }

            for (int m = table.Count; m <= n; m++) {
                v *= 2 * m - 1;
                table.Add(v);
            }

            return table[n];
        }
    }
}
