using MultiPrecision;
using System;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            const double threshold = 10.75d;

            for (MultiPrecision<Pow2.N8> x = 0; x <= threshold; x += 1 / 64d) {
                MultiPrecision<Pow2.N8> y = NearZero<Pow2.N8, Pow2.N16>.FresnelC(x);

                Console.WriteLine($"{x},{y}");
            }

            for (MultiPrecision<Pow2.N8> x = threshold; x <= 128; x += 1 / 16d) {
                MultiPrecision<Pow2.N8> y = Limit<Pow2.N8>.Fresnel(x).c;

                Console.WriteLine($"{x},{y}");
            }

            for (MultiPrecision<Pow2.N8> x = 0; x <= threshold; x += 1 / 64d) {
                MultiPrecision<Pow2.N8> y = NearZero<Pow2.N8, Pow2.N16>.FresnelS(x);

                Console.WriteLine($"{x},{y}");
            }

            for (MultiPrecision<Pow2.N8> x = threshold; x <= 128; x += 1 / 16d) {
                MultiPrecision<Pow2.N8> y = Limit<Pow2.N8>.Fresnel(x).s;

                Console.WriteLine($"{x},{y}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
