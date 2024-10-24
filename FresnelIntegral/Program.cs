using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace FresnelIntegral {
    internal class Program {
        static void Main() {
            static (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) invg16(MultiPrecision<Pow2.N16> x) {
                (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) = FGExpandN16.Value(x);

                f *= x;
                g *= x * x * x;

                return (f, g);
            };

            static (MultiPrecision<Pow2.N32> f, MultiPrecision<Pow2.N32> g) invg32(MultiPrecision<Pow2.N32> x) {
                (MultiPrecision<Pow2.N32> f, MultiPrecision<Pow2.N32> g) = FGExpandN32.Value(x);

                f *= x;
                g *= x * x * x;

                return (f, g);
            };

            using (StreamWriter sw = new("../../../../results/fresnel_scaled_fg.csv")) {
                sw.WriteLine("x,f,g");

                for (double x = 0; x <= 16; x += 1d / 32) {
                    (MultiPrecision<Pow2.N16> f, MultiPrecision<Pow2.N16> g) = invg16(x);

                    sw.WriteLine($"{x},{f},{g}");

                    Console.WriteLine($"{x}\t{f}\t{g}");
                }
            }

            List<(MultiPrecision<Pow2.N64> xmin, MultiPrecision<Pow2.N64> xmax, MultiPrecision<Pow2.N64> limit_range)> ranges = [
                (0.5, 1, 1 / 4096d), (1, 2, 1 / 4096d), (2, 4, 1 / 4096d), (4, 8, 1 / 4096d), (8, 16, 1 / 4096d), (16, 32, 1 / 4096d)
            ];

            using StreamWriter f_sw = new("../../../../results_disused/fresnel_scaled_f_pade_table.csv");
            using StreamWriter g_sw = new("../../../../results_disused/fresnel_scaled_g_pade_table.csv");

            bool approximate(MultiPrecision<Pow2.N64> xmin, MultiPrecision<Pow2.N64> xmax) {
                Console.WriteLine($"[{xmin}, {xmax}]");

                List<(MultiPrecision<Pow2.N64> x, MultiPrecision<Pow2.N64> f, MultiPrecision<Pow2.N64> g)> expecteds_range = [];

                for (MultiPrecision<Pow2.N64> x = xmin, h = (xmax - xmin) / 8192; x <= xmax; x += h) {
                    (MultiPrecision<Pow2.N32> f, MultiPrecision<Pow2.N32> g) = invg32(x.Convert<Pow2.N32>());

                    expecteds_range.Add((x, f.Convert<Pow2.N64>(), g.Convert<Pow2.N64>()));
                }

                Console.WriteLine("expecteds computed");

                MultiPrecision<Pow2.N64> f0 = expecteds_range.First().f;
                MultiPrecision<Pow2.N64> g0 = expecteds_range.First().g;

                Vector<Pow2.N64> xs = expecteds_range.Select(item => item.x - xmin).ToArray();
                Vector<Pow2.N64> fs = expecteds_range.Select(item => item.f).ToArray();
                Vector<Pow2.N64> gs = expecteds_range.Select(item => item.g).ToArray();

                SumTable<Pow2.N64> f_sum_table = new(xs, fs), g_sum_table = new(xs, gs);

                for (int coefs = 5; coefs <= 72; coefs++) {
                    foreach ((int m, int n) in CurveFittingUtils.EnumeratePadeDegree(coefs, 0)) {
                        PadeFitter<Pow2.N64> f_pade = new(f_sum_table, m, n);
                        PadeFitter<Pow2.N64> g_pade = new(g_sum_table, m, n);

                        Vector<Pow2.N64> f_param = f_pade.Fit();
                        Vector<Pow2.N64> f_errs = f_pade.Error(f_param);

                        Vector<Pow2.N64> g_param = g_pade.Fit();
                        Vector<Pow2.N64> g_errs = g_pade.Error(g_param);

                        MultiPrecision<Pow2.N64> max_rateerr = MultiPrecision<Pow2.N64>.Max(
                            CurveFittingUtils.MaxRelativeError(fs, f_pade.Regress(xs, f_param)),
                            CurveFittingUtils.MaxRelativeError(gs, g_pade.Regress(xs, g_param))
                        );

                        Console.WriteLine($"m={m},n={n}");
                        Console.WriteLine($"{max_rateerr:e20}");

                        if (max_rateerr > "1e-22") {
                            coefs += 4;
                            break;
                        }

                        if (max_rateerr < "1e-34") {
                            return false;
                        }

                        if (max_rateerr < "1e-31" &&
                            !CurveFittingUtils.HasLossDigitsPolynomialCoef(f_param[..m], 0, xmax - xmin) &&
                            !CurveFittingUtils.HasLossDigitsPolynomialCoef(f_param[m..], 0, xmax - xmin) &&
                            !CurveFittingUtils.HasLossDigitsPolynomialCoef(g_param[..m], 0, xmax - xmin) &&
                            !CurveFittingUtils.HasLossDigitsPolynomialCoef(g_param[m..], 0, xmax - xmin)) {

                            f_sw.WriteLine($"x=[{xmin},{xmax}]");
                            f_sw.WriteLine($"samples={expecteds_range.Count}");
                            f_sw.WriteLine($"m={m},n={n}");
                            f_sw.WriteLine("numer");
                            foreach (var (_, val) in f_param[..m]) {
                                f_sw.WriteLine($"{val:e38}");
                            }
                            f_sw.WriteLine("denom");
                            foreach (var (_, val) in f_param[m..]) {
                                f_sw.WriteLine($"{val:e38}");
                            }

                            f_sw.WriteLine("coef");
                            foreach ((var numer, var denom) in CurveFittingUtils.EnumeratePadeCoef(f_param, m, n)) {
                                f_sw.WriteLine($"({ToFP128(numer)}, {ToFP128(denom)}),");
                            }

                            f_sw.WriteLine("relative err");
                            f_sw.WriteLine($"{max_rateerr:e20}");
                            f_sw.Flush();

                            g_sw.WriteLine($"x=[{xmin},{xmax}]");
                            g_sw.WriteLine($"samples={expecteds_range.Count}");
                            g_sw.WriteLine($"m={m},n={n}");
                            g_sw.WriteLine("numer");
                            foreach (var (_, val) in g_param[..m]) {
                                g_sw.WriteLine($"{val:e38}");
                            }
                            g_sw.WriteLine("denom");
                            foreach (var (_, val) in g_param[m..]) {
                                g_sw.WriteLine($"{val:e38}");
                            }

                            g_sw.WriteLine("coef");
                            foreach ((var numer, var denom) in CurveFittingUtils.EnumeratePadeCoef(g_param, m, n)) {
                                g_sw.WriteLine($"({ToFP128(numer)}, {ToFP128(denom)}),");
                            }

                            g_sw.WriteLine("relative err");
                            g_sw.WriteLine($"{max_rateerr:e20}");
                            g_sw.Flush();

                            return true;
                        }
                    }
                }

                return false;
            }

            Segmenter<Pow2.N64> segmenter = new(ranges, approximate);
            segmenter.Execute();

            foreach ((var xmin, var xmax, bool is_successs) in segmenter.ApproximatedRanges) {
                f_sw.WriteLine($"[{xmin},{xmax}],{(is_successs ? "OK" : "NG")}");
            }

            foreach ((var xmin, var xmax, bool is_successs) in segmenter.ApproximatedRanges) {
                g_sw.WriteLine($"[{xmin},{xmax}],{(is_successs ? "OK" : "NG")}");
            }

            Console.WriteLine("END");
            Console.Read();
        }

        public static string ToFP128(MultiPrecision<Pow2.N64> x) {
            Sign sign = x.Sign;
            long exponent = x.Exponent;
            uint[] mantissa = x.Mantissa.Reverse().ToArray();

            string code = $"({(sign == Sign.Plus ? "+1" : "-1")}, {exponent}, 0x{mantissa[0]:X8}{mantissa[1]:X8}uL, 0x{mantissa[2]:X8}{mantissa[3]:X8}uL)";

            return code;
        }
    }
}
