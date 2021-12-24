using MultiPrecision;
using System;
using System.Collections.Generic;
using System.IO;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            using (StreamWriter sw = new("../../../../results_disused/fresnel_fg16.csv")) {
                sw.WriteLine("x,f(x),g(x)");
                for (MultiPrecision<Pow2.N16> x = 0; x <= 128; x += 1d / 64) {
                    (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) = FGExpandN16.Value(x);

                    sw.WriteLine($"{x},{f},{g}");
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
