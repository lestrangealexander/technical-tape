/* Virtual tape: the per-sample core (runs at the oversampled rate).
   Memory nonlinearity: tanh with a low-passed feedback of its own output (level-dependent phase and
   frequency response), a small even-order term (asymmetry), then dynamic self-erasure: a one-pole
   low-pass whose cutoff falls as the high-frequency envelope of the record signal rises. */
#include <math.h>

void vt_process(const double* u, const double* uhf, double* out, long n,
                double g, double kappa, double fbCoef, double a2,
                double f0, double e0, double envRel, double fsOs)
{
    double s = 0.0, env = 0.0, y = 0.0;
    const double twoPiOverFs = 6.283185307179586 / fsOs;
    for (long i = 0; i < n; ++i)
    {
        const double v = g * (u[i] - kappa * s);
        const double t = tanh (v) / g;
        const double m = t + a2 * g * t * t;
        s += fbCoef * (m - s);

        const double h = fabs (uhf[i]);
        env = h > env ? h : env * envRel;
        const double r = env / e0;
        const double fc = f0 / (1.0 + r * r);
        const double a = 1.0 - exp (-twoPiOverFs * fc);
        y += a * (m - y);
        out[i] = y;
    }
}
