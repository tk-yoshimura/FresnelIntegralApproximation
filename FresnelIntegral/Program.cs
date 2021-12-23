using MultiPrecision;
using System;
using System.Collections.Generic;
using System.IO;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            using (StreamWriter sw = new("../../../../results/fresnelc_point5s.csv")) {
                sw.WriteLine("index,x");

                List<MultiPrecision<Pow2.N8>> xs = new() {
                    FindPoint5N8.FresnelC(0.5),
                    FindPoint5N8.FresnelC(1.4),
                    FindPoint5N8.FresnelC(2.0),
                    FindPoint5N8.FresnelC(2.4),
                };

                for (int i = 0; i < xs.Count; i++) {
                    sw.WriteLine($"{i},{xs[i]}");
                }

                while (xs.Count < 4096) {
                    MultiPrecision<Pow2.N8> x0 = xs[^3], x1 = xs[^2], x2 = xs[^1];
                    MultiPrecision<Pow2.N8> x_predict = x2 + MultiPrecision<Pow2.N8>.Square(x2 - x1) / (x1 - x0);

                    MultiPrecision<Pow2.N8> x = FindPoint5N8.FresnelC(x_predict);
                    MultiPrecision<Pow2.N8> y = FresnelN8.FresnelC(x);
                    MultiPrecision<Pow2.N8> y_dec = FresnelN8.FresnelC(x - 1e-74) - 0.5d;
                    MultiPrecision<Pow2.N8> y_inc = FresnelN8.FresnelC(x + 1e-74) - 0.5d;

                    sw.WriteLine($"{xs.Count},{x}");
                    xs.Add(x);

                    Console.WriteLine($"{x}\t{y}\t{y_dec:e3}\t{y_inc:e3}");
                }
            }

            using (StreamWriter sw = new("../../../../results/fresnels_point5s.csv")) {
                sw.WriteLine("index,x");

                List<MultiPrecision<Pow2.N8>> xs = new() {
                    FindPoint5N8.FresnelS(1.1),
                    FindPoint5N8.FresnelS(1.7),
                    FindPoint5N8.FresnelS(2.2),
                    FindPoint5N8.FresnelS(2.6),
                };

                for (int i = 0; i < xs.Count; i++) {
                    sw.WriteLine($"{i},{xs[i]}");
                }

                while (xs.Count < 4096) {
                    MultiPrecision<Pow2.N8> x0 = xs[^3], x1 = xs[^2], x2 = xs[^1];
                    MultiPrecision<Pow2.N8> x_predict = x2 + MultiPrecision<Pow2.N8>.Square(x2 - x1) / (x1 - x0);

                    MultiPrecision<Pow2.N8> x = FindPoint5N8.FresnelS(x_predict);
                    MultiPrecision<Pow2.N8> y = FresnelN8.FresnelS(x);
                    MultiPrecision<Pow2.N8> y_dec = FresnelN8.FresnelS(x - 1e-74) - 0.5d;
                    MultiPrecision<Pow2.N8> y_inc = FresnelN8.FresnelS(x + 1e-74) - 0.5d;

                    sw.WriteLine($"{xs.Count},{x}");
                    xs.Add(x);

                    Console.WriteLine($"{x}\t{y}\t{y_dec:e3}\t{y_inc:e3}");
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
