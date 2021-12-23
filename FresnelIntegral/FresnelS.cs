using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MultiPrecision;

namespace FresnelIntegral {
    internal static class FresnelS<N, M> where N : struct, IConstant where M: struct, IConstant {
        public static MultiPrecision<N> Value(MultiPrecision<N> x) {
            if (x.Sign == Sign.Minus) {
                return -Value(-x);
            }

            MultiPrecision<N> v = MultiPrecision<N>.PI
        }
    }
}
