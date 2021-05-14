from operator import mul, pow, sub, truediv

from util.func_utils import nest, part, pipe, reduce


# FUNCTIONS


def inv(x):
    return 1 / x


def isnan(x):
    return x != x


def abserr(a, b):
    return abs(a - b)


def relerr(a, b):
    return abserr(a, b) / a


def mean(itr):
    n = len(itr)
    if (n <= 0):
        return float("nan")
    elif (n == 1):
        return next(iter(itr))
    return sum(itr) / len(itr)


def harmonic_mean(itr):
    return inv(mean(pipe(inv)(itr)))


def geometric_mean(itr):
    n = len(itr)
    if (n <= 0):
        return float("nan")
    elif (n == 1):
        return next(iter(itr))
    return reduce(mul)(itr) ** inv(len(itr))


def median(itr):
    n = len(itr)
    m = n // 2
    srt = sorted(itr)
    return srt[m] if n % 2 else (srt[m - 1] + srt[m]) / 2


def mode(itr):
    dct = dict()
    for n in itr:
        dct[n] = dct[n] + 1 if (n in dct) else 1
    return max(dct, key=dct.get)


def multimode(itr):
    dct = dict()
    for n in itr:
        dct[n] = dct[n] + 1 if (n in dct) else 1
    x = max(dct.values())
    return {k for k, v in dct.items() if v == x}


def var(itr):
    if len(itr) <= 1:
        return float("nan")
    mu = mean(itr)
    op = nest(part(sub, mu), part(pow, 2))
    return sum(pipe(op)(itr)) / (len(itr) - 1)


def sd(itr):
    return var(itr) ** 0.5


def zscore(itr):
    mu = mean(itr)
    s = sd(itr)
    op = nest(part(sub, mu), part(truediv, s))
    return pipe(op)(itr)


def minmax(itr):
    if len(itr) <= 1:
        return
    lo, hi = min(itr), max(itr)
    rng = (hi - lo)
    op = nest(part(sub, lo), part(truediv, rng))
    return pipe(op)(itr)
