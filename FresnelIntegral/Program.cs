using MultiPrecision;
using System;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            MultiPrecision<Pow2.N16> x = 2;

            (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) = FGExpandN16.Value(x);

            (MultiPrecision<Pow2.N16>[] fs, MultiPrecision<Pow2.N16>[] gs) = FGDiffSeries<Pow2.N16>.Diff(x, f, g, 6);

            for (int i = 0; i < fs.Length; i++) {
                Console.WriteLine(fs[i]);
            }

            for (int i = 0; i < gs.Length; i++) {
                Console.WriteLine(gs[i]);
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
