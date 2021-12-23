using MultiPrecision;
using System;

namespace FresnelIntegral {
    internal class Program {
        static void Main(string[] args) {
            for (MultiPrecision<Pow2.N32> x = 0; x <= 16; x += 1 / 64d) {
                MultiPrecision<Pow2.N32> y = FresnelC<Pow2.N32, Pow2.N64>.Value(x);

                Console.WriteLine($"{x},{y}");

                if (!(y < 1)) {
                    break;
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
