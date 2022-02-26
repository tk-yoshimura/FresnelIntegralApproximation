using MultiPrecision;
using System;
using System.Collections.Generic;
using System.IO;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            for (double x = 16; x <= 20; x += 1d / 16) {
                MultiPrecision<Pow2.N16> y = FresnelN16.FresnelS(x);

                Console.WriteLine($"\"{y:e40}\",");
            }


            static (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) invg16(MultiPrecision<Pow2.N16> x) {
                MultiPrecision<Pow2.N16> v = MultiPrecision<Pow2.N16>.Pow2(x);
                (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) = FGExpandN16.Value(v);

                f = MultiPrecision<Pow2.N16>.Log2(f);
                g = MultiPrecision<Pow2.N16>.Log2(g);
            
                return (f, g);
            };

             static (MultiPrecision<Pow2.N32> f, MultiPrecision<Pow2.N32> g) invg32(MultiPrecision<Pow2.N32> x) {
                MultiPrecision<Pow2.N32> v = MultiPrecision<Pow2.N32>.Pow2(x);
                (MultiPrecision<Pow2.N32> f, MultiPrecision<Pow2.N32> g) = FGExpandN32.Value(v);

                f = MultiPrecision<Pow2.N32>.Log2(f);
                g = MultiPrecision<Pow2.N32>.Log2(g);
            
                return (f, g);
            };

            using (StreamWriter sw = new("../../../../results/fresnel_log2log2_fg.csv")) {
                sw.WriteLine("x,f,g");

                for (double x = 0; x <= 8; x += 1d / 32) {
                    (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) = invg16(x);

                    sw.WriteLine($"{x},{f},{g}");

                    Console.WriteLine($"{x}\t{f}\t{g}");
                }
            }
                        
            using (StreamWriter sw = new("../../../../results/fresnel_log2log2_pade_table.csv")) {
                for (MultiPrecision<Pow2.N16> x = 0; x <= 4; x += 0.5d) {

                    bool is_e31 = false;

                    MultiPrecision<Pow2.N16> ddx = Math.ScaleB(1, -2);
                    MultiPrecision<Pow2.N32>[] fds = FiniteDifference<Pow2.N32>.Diff(
                        x.Convert<Pow2.N32>(), (x) => invg32(x).f, Math.ScaleB(1, -24)
                    );
                    MultiPrecision<Pow2.N32>[] gds = FiniteDifference<Pow2.N32>.Diff(
                        x.Convert<Pow2.N32>(), (x) => invg32(x).g, Math.ScaleB(1, -24)
                    );

                    List<MultiPrecision<Pow2.N16>> fexpecteds = new(), gexpecteds = new();

                    for (MultiPrecision<Pow2.N16> dx = -ddx, h = ddx / 512; dx <= ddx; dx += h) {
                        (MultiPrecision<Pow2.N16> fexpected, MultiPrecision<Pow2.N16> gexpected)
                            = invg16(x + dx);

                        fexpecteds.Add(fexpected);
                        gexpecteds.Add(gexpected);
                    }

                    for (int n = 4; n <= 16; n += 1) {
                        MultiPrecision<Pow2.N16>[] fcs = new MultiPrecision<Pow2.N16>[n * 2 + 1];
                        MultiPrecision<Pow2.N16>[] gcs = new MultiPrecision<Pow2.N16>[n * 2 + 1];
                            
                        (fcs[0], gcs[0]) = invg16(x);
                        for (int i = 0; i < n * 2; i++) {
                            fcs[i + 1] = fds[i].Convert<Pow2.N16>() * MultiPrecision<Pow2.N16>.TaylorSequence[i + 1];
                        }
                        for (int i = 0; i < n * 2; i++) {
                            gcs[i + 1] = gds[i].Convert<Pow2.N16>() * MultiPrecision<Pow2.N16>.TaylorSequence[i + 1];
                        }

                        (MultiPrecision<Pow2.N16>[] fms, MultiPrecision<Pow2.N16>[] fns) = PadeSolver<Pow2.N16>.Solve(fcs, n, n);
                        (MultiPrecision<Pow2.N16>[] gms, MultiPrecision<Pow2.N16>[] gns) = PadeSolver<Pow2.N16>.Solve(gcs, n, n);

                        MultiPrecision<Pow2.N16> err = 0;

                        for ((MultiPrecision<Pow2.N16> dx, MultiPrecision<Pow2.N16> h, int i) = (-ddx, ddx / 512, 0); i < fexpecteds.Count; dx += h, i++) {
                            MultiPrecision<Pow2.N16> expected = fexpecteds[i];
                            MultiPrecision<Pow2.N16> actual = PadeSolver<Pow2.N16>.Approx(dx, fms, fns);

                            err = MultiPrecision<Pow2.N16>.Max(err, MultiPrecision<Pow2.N16>.Abs(expected / actual - 1));
                        }

                        for ((MultiPrecision<Pow2.N16> dx, MultiPrecision<Pow2.N16> h, int i) = (-ddx, ddx / 512, 0); i < gexpecteds.Count; dx += h, i++) {
                            MultiPrecision<Pow2.N16> expected = gexpecteds[i];
                            MultiPrecision<Pow2.N16> actual = PadeSolver<Pow2.N16>.Approx(dx, gms, gns);

                            err = MultiPrecision<Pow2.N16>.Max(err, MultiPrecision<Pow2.N16>.Abs(expected / actual - 1));
                        }

                        Console.WriteLine($"x={x}, n={n}, |dx| = {ddx}");
                        Console.WriteLine($"relative error = {err:e10}");

                        if (err < 2e-31) {
                            sw.WriteLine($"x={x}, n={n}, |dx| = {ddx}");

                            sw.WriteLine($",i,p_i,q_i");
                            for (int i = 0; i <= n; i++) {
                                sw.WriteLine($"f,{i},{fms[i]:e64},{fns[i]:e64}");
                            }
                            for (int i = 0; i <= n; i++) {
                                sw.WriteLine($"g,{i},{gms[i]:e64},{gns[i]:e64}");
                            }

                            sw.WriteLine($"relative error = {err:e10}\n");

                            is_e31 = true;

                            sw.Flush();

                            break;
                        }
                    }

                    if (!is_e31) {
                        sw.WriteLine($"{x} not convergence.");
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
