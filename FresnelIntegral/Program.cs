using MultiPrecision;
using System;
using System.Collections.Generic;
using System.IO;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            using (StreamWriter sw = new("../../../../results/fresnel_fg.csv")) {
                sw.WriteLine("x,f(x),g(x)");
                for (MultiPrecision<Pow2.N8> x = 0; x <= 128; x += 1d / 64) {
                    (MultiPrecision<Pow2.N8> f, MultiPrecision<Pow2.N8> g) = FGExpandN8.Value(x);

                    sw.WriteLine($"{x},{f},{g}");
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
