using MultiPrecision;
using System;
using System.IO;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            MultiPrecision<Pow2.N16> bin = MultiPrecision<Pow2.N16>.Ldexp(1, -2);

            using (StreamWriter sw = new("../../../../results/pade_table.csv")) {

                for (MultiPrecision<Pow2.N16> x = bin; x <= 12; x += bin) {
                    (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) = FGExpandN16.Value(x);

                    for (int n = 4; n <= 24; n += 1) {
                        sw.WriteLine($"x={x}, n={n}");
                        sw.WriteLine($",i,p_i,q_i");

                        (MultiPrecision<Pow2.N16>[] fs, MultiPrecision<Pow2.N16>[] gs) = FGDiffSeries<Pow2.N16>.Diff(x, f, g, 2 * n);

                        for (int i = 0; i <= 2 * n; i++) {
                            fs[i] *= MultiPrecision<Pow2.N16>.TaylorSequence[i];
                            gs[i] *= MultiPrecision<Pow2.N16>.TaylorSequence[i];
                        }

                        (MultiPrecision<Pow2.N16>[] fms, MultiPrecision<Pow2.N16>[] fns) = PadeSolver<Pow2.N16>.Solve(fs, n, n);
                        (MultiPrecision<Pow2.N16>[] gms, MultiPrecision<Pow2.N16>[] gns) = PadeSolver<Pow2.N16>.Solve(gs, n, n);

                        for (int i = 0; i <= n; i++) {
                            sw.WriteLine($"f(x),{i},{fms[i]:e64},{fns[i]:e64}");
                        }
                        for (int i = 0; i <= n; i++) {
                            sw.WriteLine($"g(x),{i},{gms[i]:e64},{gns[i]:e64}");
                        }

                        MultiPrecision<Pow2.N16> f_err = 0, g_err = 0;

                        for (MultiPrecision<Pow2.N16> dx = -bin / 2, ddx = bin / 512; dx <= bin / 2; dx += ddx) {
                            if (x + dx <= 0) {
                                continue;
                            }

                            (MultiPrecision<Pow2.N16> f_expected, MultiPrecision<Pow2.N16> g_expected) = FGExpandN16.Value(x + dx);

                            MultiPrecision<Pow2.N16> f_actual = PadeSolver<Pow2.N16>.Approx(dx, fms, fns);
                            MultiPrecision<Pow2.N16> g_actual = PadeSolver<Pow2.N16>.Approx(dx, gms, gns);

                            f_err = MultiPrecision<Pow2.N16>.Max(f_err, MultiPrecision<Pow2.N16>.Abs(f_expected / f_actual - 1));
                            g_err = MultiPrecision<Pow2.N16>.Max(g_err, MultiPrecision<Pow2.N16>.Abs(g_expected / g_actual - 1));
                        }

                        sw.WriteLine($"relative error f(x)={f_err:e10}, relative error g(x)={g_err:e10}\n");

                        Console.WriteLine($"x={x}, n={n}");
                        Console.WriteLine($"relative error f(x)={f_err:e10}");
                        Console.WriteLine($"relative error g(x)={g_err:e10}");

                        if (f_err < 1e-34 && g_err < 1e-34) {
                            break;
                        }
                    }
                }
            }

            //using (StreamWriter sw = new("../../../../results_disused/fresnelc_table.csv")) {
            //    for (double x = 1d / 16; x <= 16; x += 1d / 16) {
            //        sw.WriteLine($"\"{FresnelN8.FresnelC(x):e40}\",");
            //    }
            //}
            //
            //using (StreamWriter sw = new("../../../../results_disused/fresnels_table.csv")) {
            //    for (double x = 1d / 16; x <= 16; x += 1d / 16) {
            //        sw.WriteLine($"\"{FresnelN8.FresnelS(x):e40}\",");
            //    }
            //}

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
