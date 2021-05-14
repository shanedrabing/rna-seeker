import functools
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor


# FUNCTIONS


def identity(x):
    return x


def null(x):
    pass


def transpose(x):
    return tuple(zip(*x))


def getitems(index, itr):
    return tuple(map(itr.__getitem__, index))


# FUNCTIONS (WRAPPERS)


def memoize(f):
    memo = dict()

    @functools.wraps(f)
    def memoizef(*args, **kwargs):
        key = "%s%s" % (args, kwargs)
        if key not in memo:
            memo[key] = f(*args, **kwargs)
        return memo[key]
    return memoizef


# FUNCTIONS (FUNCTIONAL-PROGRAMMING)


def accum(f):
    def accumf(itr):
        x, *itr = itr
        yield x
        for y in itr:
            x = f(x, y)
            yield x
    return accumf


def apply(f):
    def applyf(*itrs):
        return map(f, *itrs)
    return applyf


def call(x):
    def callf(f):
        return f(x)
    return callf


def catch(f, err, cf=null):
    def catchf(x):
        try:
            return f(x)
        except err:
            return cf(x)
    return catchf


def chunker(n):
    def chunkerf(itr):
        while itr:
            x, itr = itr[:n], itr[n:]
            yield x
    return chunkerf


def display(f=identity):
    def displayf(x):
        print(f(x))
        return x
    return displayf


def filt(f):
    def filtf(x):
        return filter(f, x)
    return filtf


def get(x, iskey=True):
    def getf(y):
        try:
            return y[x] if iskey else x[y]
        except KeyError:
            return None
    return getf


def log(message, *attrs):
    log.count = 0
    def logf(*_):
        log.count += 1
        return message.format(*map(getattr, (log,) * len(attrs), attrs))
    return logf


def nest(*fs):
    def nestf(*x):
        for i, f in enumerate(fs):
            x = f(*x) if (i == 0) else f(x)
        return x
    return nestf


def pack(f):
    def packf(*args):
        return f(args)
    return packf


def part(f, *args, **kwargs):
    def partf(*x):
        return f(*x, *args, **kwargs)
    return partf


def pipe(f):
    def pipef(*x):
        return tuple(map(f, *x))
    return pipef


def reduce(f):
    def reducef(itr):
        x, *itr = itr
        for y in itr:
            x = f(x, y)
        return x
    return reducef


def unpack(f):
    def unpackf(itr):
        return f(*itr)
    return unpackf


# FUNCTIONS (MULTI-THREADING)


def ripper(f, exectype="thread"):
    if (exectype == "thread"):
        exectype = ThreadPoolExecutor
    elif (exectype == "process"):
        exectype = ProcessPoolExecutor

    def ripperf(*x):
        with exectype() as executor:
            return executor.map(f, *x)
    return ripperf
