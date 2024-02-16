from functools import reduce
from pathlib import Path
from typing import Sequence

from py_ecc import optimized_bn128 as bn128

from snarkpy.field import GF
from snarkpy.hash import keccak
from snarkpy.plonk.parse import parse_vkey
from snarkpy.plonk.prove import Proof


def verify(
    vkey_path: Path,
    public_inputs_int: Sequence[int],
    proof: Proof,
    debug: bool = False,
) -> bool:
    vkey = parse_vkey(vkey_path)
    assert vkey.curve == "bn128"
    Fr = GF(bn128.curve_order, 2**256)
    omega = Fr(vkey.omega)
    assert vkey.n_public == len(public_inputs_int)
    public_inputs = [Fr(pi) for pi in public_inputs_int]
    transcript1 = b""
    transcript1 += int(bn128.normalize(vkey.qm)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.qm)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.ql)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.ql)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.qr)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.qr)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.qo)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.qo)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.qc)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.qc)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.s1)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.s1)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.s2)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.s2)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.s3)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(vkey.s3)[1]).to_bytes(32, "big")
    for i in range(vkey.n_public):
        transcript1 += Fr.to_bytes(public_inputs[i])[::-1]
    transcript1 += int(bn128.normalize(proof.a)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(proof.a)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(proof.b)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(proof.b)[1]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(proof.c)[0]).to_bytes(32, "big")
    transcript1 += int(bn128.normalize(proof.c)[1]).to_bytes(32, "big")
    h = keccak(transcript1)
    beta = Fr(int.from_bytes(h, "big"))
    transcript2 = Fr.to_bytes(beta)[::-1]
    h = keccak(transcript2)
    gamma = Fr(int.from_bytes(h, "big"))

    transcript3 = b""
    transcript3 += Fr.to_bytes(beta)[::-1]
    transcript3 += Fr.to_bytes(gamma)[::-1]
    transcript3 += int(bn128.normalize(proof.z)[0]).to_bytes(32, "big")
    transcript3 += int(bn128.normalize(proof.z)[1]).to_bytes(32, "big")
    h = keccak(transcript3)
    alpha = Fr(int.from_bytes(h, "big"))

    transcript4 = b""
    transcript4 += Fr.to_bytes(alpha)[::-1]
    transcript4 += int(bn128.normalize(proof.t_low)[0]).to_bytes(32, "big")
    transcript4 += int(bn128.normalize(proof.t_low)[1]).to_bytes(32, "big")
    transcript4 += int(bn128.normalize(proof.t_mid)[0]).to_bytes(32, "big")
    transcript4 += int(bn128.normalize(proof.t_mid)[1]).to_bytes(32, "big")
    transcript4 += int(bn128.normalize(proof.t_high)[0]).to_bytes(32, "big")
    transcript4 += int(bn128.normalize(proof.t_high)[1]).to_bytes(32, "big")
    h = keccak(transcript4)
    zeta = Fr(int.from_bytes(h, "big"))

    transcript5 = b""
    transcript5 += Fr.to_bytes(zeta)[::-1]
    transcript5 += int(proof.eval_a).to_bytes(32, "big")
    transcript5 += int(proof.eval_b).to_bytes(32, "big")
    transcript5 += int(proof.eval_c).to_bytes(32, "big")
    transcript5 += int(proof.eval_s1).to_bytes(32, "big")
    transcript5 += int(proof.eval_s2).to_bytes(32, "big")
    transcript5 += int(proof.eval_zw).to_bytes(32, "big")
    h = keccak(transcript5)
    v = Fr(int.from_bytes(h, "big"))

    transcript6 = b""
    transcript6 += Fr.to_bytes(v)[::-1]
    transcript6 += int(bn128.normalize(proof.w_zeta)[0]).to_bytes(32, "big")
    transcript6 += int(bn128.normalize(proof.w_zeta)[1]).to_bytes(32, "big")
    transcript6 += int(bn128.normalize(proof.w_zetaw)[0]).to_bytes(32, "big")
    transcript6 += int(bn128.normalize(proof.w_zetaw)[1]).to_bytes(32, "big")
    h = keccak(transcript6)
    u = Fr(int.from_bytes(h, "big"))

    if debug:
        print(f"{public_inputs = }")
        print(f"{beta = }")
        print(f"{gamma = }")
        print(f"{alpha = }")
        print(f"{zeta = }")
        print(f"{v = }")
        print(f"{u = }")

    n = 2**vkey.power
    eval_zh = zeta**n - Fr(1)
    eval_ls = []
    for i in range(max(vkey.n_public, 1)):
        eval_ls.append(omega**i * (zeta**n - Fr(1)) / Fr(n) / (zeta - omega**i))
    eval_l1 = (zeta**n - Fr(1)) / Fr(n) / (zeta - Fr(1))
    eval_pi = sum(
        [public_inputs[i] * eval_ls[i] for i in range(vkey.n_public)], Fr(0)
    ) * Fr(-1)

    eval_a = proof.eval_a
    eval_b = proof.eval_b
    eval_c = proof.eval_c
    eval_s1 = proof.eval_s1
    eval_s2 = proof.eval_s2
    eval_zw = proof.eval_zw

    r0 = (
        eval_pi
        - eval_l1 * alpha**2
        - alpha
        * (eval_a + beta * eval_s1 + gamma)
        * (eval_b + beta * eval_s2 + gamma)
        * (eval_c + gamma)
        * eval_zw
    )

    d = reduce(
        lambda x, y: bn128.add(x, y),
        [
            bn128.multiply(vkey.qm, int(eval_a * eval_b)),
            bn128.multiply(vkey.ql, int(eval_a)),
            bn128.multiply(vkey.qr, int(eval_b)),
            bn128.multiply(vkey.qo, int(eval_c)),
            vkey.qc,
            bn128.multiply(
                proof.z,
                int(
                    (eval_a + beta * zeta + gamma)
                    * (eval_b + beta * vkey.k1 * zeta + gamma)
                    * (eval_c + beta * vkey.k2 * zeta + gamma)
                    * alpha
                    + eval_l1 * alpha**2
                    + u
                ),
            ),
            bn128.multiply(
                vkey.s3,
                int(
                    (eval_a + beta * eval_s1 + gamma)
                    * (eval_b + beta * eval_s2 + gamma)
                    * alpha
                    * beta
                    * eval_zw
                    * Fr(-1)
                ),
            ),
            bn128.multiply(proof.t_low, int(eval_zh * Fr(-1))),
            bn128.multiply(proof.t_mid, int(zeta**n * eval_zh * Fr(-1))),
            bn128.multiply(proof.t_high, int(zeta ** (2 * n) * eval_zh * Fr(-1))),
        ],
    )
    f = reduce(
        lambda x, y: bn128.add(x, y),
        [
            d,
            bn128.multiply(proof.a, int(v)),
            bn128.multiply(proof.b, int(v**2)),
            bn128.multiply(proof.c, int(v**3)),
            bn128.multiply(vkey.s1, int(v**4)),
            bn128.multiply(vkey.s2, int(v**5)),
        ],
    )
    e = bn128.multiply(
        bn128.G1,
        int(
            r0 * Fr(-1)
            + v * eval_a
            + v**2 * eval_b
            + v**3 * eval_c
            + v**4 * eval_s1
            + v**5 * eval_s2
            + u * eval_zw
        ),
    )
    e1 = bn128.pairing(
        vkey.x_2, bn128.add(proof.w_zeta, bn128.multiply(proof.w_zetaw, int(u)))
    )
    e2 = bn128.pairing(
        bn128.G2,
        reduce(
            lambda x, y: bn128.add(x, y),
            [
                bn128.multiply(proof.w_zeta, int(zeta)),
                bn128.multiply(proof.w_zetaw, int(u * zeta * omega)),
                f,
                bn128.neg(e),
            ],
        ),
    )
    return e1 == e2  # verification
