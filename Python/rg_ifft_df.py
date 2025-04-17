import numpy as np
from scipy import interpolate
import scipy.optimize as opt


def f1(r0, param):
    ss = 45 / 0.3
    ro = 10 ** r0
    kappa = param[3] / (ro[1] - ro[0]) ** 2
    jj = 2 * kappa
    lfr = np.zeros(len(ro))
    lrr = np.zeros(len(ro))
    for i in range(len(ro)):
        ll = (ro[1] - ro[0]) * param[i + 4]
        r1 = ro[i] - ll
        r2 = ro[i] + ll

        rate = 1000
        del_r = (r2 - r1) / rate  # sampling interval

        # generate discrete points
        r = np.arange(r1, r2, del_r)
        k = np.fft.fftfreq(len(r), del_r)
        r = r[10:len(r) - 10]
        k = k[10:len(k) - 10]

        g = 1 / (4 * jj * ((k / ll ** param[2]) ** 2 / 2 - (k / ll ** param[2]) ** 4 / 24)) * \
            (np.sin(k / 2) / np.sin(k / ll ** param[2] / 2)) ** param[1]
        yy = 1 / ll ** ((2 * param[0] + 1) * param[2]) * np.fft.ifft(g).real
        for j in range(len(yy)):
            if yy[j] < 10 ** (-100):
                yy[j] = 10 ** (-100)

        lf = np.log10(yy)
        lr = np.log10(r)

        lfr[i] = lf[np.int(0.5 * len(yy))]
        lrr[i] = lr[np.int(0.5 * len(yy))]

    return lrr, lfr


def ev(ob, init, bb):
    beta = opt.least_squares(ob, init, bounds=bb)
    return beta


def th(x0, param):
    yy = np.zeros(len(x0))
    f = f1(x0, param)
    for j in range(len(x0)):
        if x0[j] > np.max(f[0]) or x0[j] < np.min(f[0]):
            fitted_curve_b = interpolate.interp1d(f[0], f[1], kind="cubic", fill_value='extrapolate')
        else:
            fitted_curve_b = interpolate.interp1d(f[0], f[1], kind="cubic")
        yy[j] = fitted_curve_b(x0[j])

    return yy
