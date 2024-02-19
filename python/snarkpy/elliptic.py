from functools import reduce
from multiprocessing import Pool
from typing import Sequence

from py_ecc import optimized_bn128 as bn128
from py_ecc.optimized_bn128 import FQ

from snarkpy.polynomial import Polynomial


def multiply(x):
    return bn128.multiply(x[0], x[1])


def exp_tau(poly: Polynomial, taus: Sequence[tuple[FQ, FQ, FQ]]) -> tuple[FQ, FQ]:
    assert poly.degree() + 1 <= len(taus)
    print(poly.degree())
    # poly.calc_coeffs_if_necessary()
    xs = [(taus[i], int(poly[i])) for i in range(poly.degree() + 1)]
    with Pool(16) as p:
        res = p.map(multiply, xs)
    return bn128.normalize(reduce(lambda x, y: bn128.add(x, y), res))
