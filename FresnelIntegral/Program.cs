using MultiPrecision;
using System;
using System.IO;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            using (StreamWriter sw = new("../../../../results/fresnel_value_table.csv")) {
                sw.WriteLine("x,fresnel_c(x),fresnel_s(x)");

                for (MultiPrecision<Pow2.N8> x = 0; x <= 128; x += 1 / 64d) {
                    (MultiPrecision<Pow2.N8> c, MultiPrecision<Pow2.N8> s) = FresnelN8.Value(x);

                    sw.WriteLine($"{x},{c},{s}");
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
