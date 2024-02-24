from contextlib import contextmanager
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from secrets import randbelow
from typing import Any, Generator, Mapping, Sequence

from py_ecc import optimized_bn128 as bn128
from py_ecc.optimized_bn128 import FQ

from snarkpy.elliptic import exp_tau
from snarkpy.field import GF, GFElement
from snarkpy.hash import keccak
from snarkpy.polynomial import Polynomial, SparsePolynomial
from snarkpy.plonk.parse import parse_wtns, parse_zkey


@contextmanager
def log_elapsed_time(
    message: str, debug: bool, indent: int = 0
) -> Generator[None, None, None]:
    prefix = " " * indent * 4
    if debug:
        print(f"{prefix}{message} start")
    now = time.time()
    yield
    elapsed = time.time() - now
    if debug:
        print(f"{prefix}{message} end ({elapsed:.2f} seconds elapsed)")


@dataclass
class Proof:
    a: tuple[FQ, FQ, FQ]
    b: tuple[FQ, FQ, FQ]
    c: tuple[FQ, FQ, FQ]
    z: tuple[FQ, FQ, FQ]
    t_low: tuple[FQ, FQ, FQ]
    t_mid: tuple[FQ, FQ, FQ]
    t_high: tuple[FQ, FQ, FQ]
    w_zeta: tuple[FQ, FQ, FQ]
    w_zetaw: tuple[FQ, FQ, FQ]
    eval_a: GFElement
    eval_b: GFElement
    eval_c: GFElement
    eval_s1: GFElement
    eval_s2: GFElement
    eval_zw: GFElement
    protocol: str = "plonk"
    curve: str = "bn128"

    def __post_init__(self):
        assert self.protocol == "plonk"
        assert self.curve == "bn128"
        Fr = GF(bn128.curve_order, 2**256)
        assert bn128.is_on_curve(self.a, bn128.b)
        assert bn128.is_on_curve(self.b, bn128.b)
        assert bn128.is_on_curve(self.c, bn128.b)
        assert bn128.is_on_curve(self.z, bn128.b)
        assert bn128.is_on_curve(self.t_low, bn128.b)
        assert bn128.is_on_curve(self.t_mid, bn128.b)
        assert bn128.is_on_curve(self.t_high, bn128.b)
        assert bn128.is_on_curve(self.w_zeta, bn128.b)
        assert bn128.is_on_curve(self.w_zetaw, bn128.b)
        assert self.eval_a in Fr
        assert self.eval_b in Fr
        assert self.eval_c in Fr
        assert self.eval_s1 in Fr
        assert self.eval_s2 in Fr
        assert self.eval_zw in Fr

    def to_json(self) -> str:
        def dict_factory(items: Sequence[tuple[str, Any]]) -> Mapping[str, Any]:
            key_map = {
                "a": "A",
                "b": "B",
                "c": "C",
                "z": "Z",
                "t_low": "T1",
                "t_mid": "T2",
                "t_high": "T3",
                "w_zeta": "Wxi",
                "w_zetaw": "Wxiw",
            }
            d = {}
            for k, v in items:
                if isinstance(v, tuple):
                    v = list(map(str, v))
                elif isinstance(v, GFElement):
                    v = str(v)
                if k in key_map:
                    k = key_map[k]
                d[k] = v
            return d

        return json.dumps(asdict(self, dict_factory=dict_factory))

    @classmethod
    def from_json(cls, res: Mapping[str, Any]) -> "Proof":
        assert res["curve"] == "bn128"
        Fr = GF(bn128.curve_order, 2**256)
        return Proof(
            a=(FQ(int(res["A"][0])), FQ(int(res["A"][1])), FQ(int(res["A"][2]))),
            b=(FQ(int(res["B"][0])), FQ(int(res["B"][1])), FQ(int(res["B"][2]))),
            c=(FQ(int(res["C"][0])), FQ(int(res["C"][1])), FQ(int(res["C"][2]))),
            z=(FQ(int(res["Z"][0])), FQ(int(res["Z"][1])), FQ(int(res["Z"][2]))),
            t_low=(FQ(int(res["T1"][0])), FQ(int(res["T1"][1])), FQ(int(res["T1"][2]))),
            t_mid=(FQ(int(res["T2"][0])), FQ(int(res["T2"][1])), FQ(int(res["T2"][2]))),
            t_high=(
                FQ(int(res["T3"][0])),
                FQ(int(res["T3"][1])),
                FQ(int(res["T3"][2])),
            ),
            w_zeta=(
                FQ(int(res["Wxi"][0])),
                FQ(int(res["Wxi"][1])),
                FQ(int(res["Wxi"][2])),
            ),
            w_zetaw=(
                FQ(int(res["Wxiw"][0])),
                FQ(int(res["Wxiw"][1])),
                FQ(int(res["Wxiw"][2])),
            ),
            eval_a=Fr(int(res["eval_a"])),
            eval_b=Fr(int(res["eval_b"])),
            eval_c=Fr(int(res["eval_c"])),
            eval_s1=Fr(int(res["eval_s1"])),
            eval_s2=Fr(int(res["eval_s2"])),
            eval_zw=Fr(int(res["eval_zw"])),
            protocol=res["protocol"],
            curve=res["curve"],
        )


def prove(
    zkey_path: Path, wtns_path: Path, debug: bool = False
) -> tuple[Proof, Sequence[int]]:
    with log_elapsed_time("parse", debug):
        zkey = parse_zkey(zkey_path)
        wtns = parse_wtns(wtns_path, zkey)

    Fr = GF(zkey.r, 2**256)

    omega = Fr.nth_root_of_unity(2**zkey.power)

    zh_coeffs = [Fr(0)] * zkey.domain_size * 4
    zh_coeffs[0] = Fr(-1)
    zh_coeffs[zkey.domain_size] = Fr(1)
    poly_zh_4 = Polynomial(Fr, coeffs=zh_coeffs)

    # round 1
    b = [Fr(randbelow(zkey.r)) for _ in range(11)]
    poly_a = wtns.poly_a
    poly_b = wtns.poly_b
    poly_c = wtns.poly_c
    taus = zkey.taus
    poly_s1 = zkey.poly_s1
    poly_s2 = zkey.poly_s2
    poly_s3 = zkey.poly_s3
    poly_l_4s = zkey.poly_l_4s
    poly_qm_4 = zkey.poly_qm_4
    poly_ql_4 = zkey.poly_ql_4
    poly_qr_4 = zkey.poly_qr_4
    poly_qo_4 = zkey.poly_qo_4
    poly_qc_4 = zkey.poly_qc_4
    poly_s1_4 = zkey.poly_s1_4
    poly_s2_4 = zkey.poly_s2_4
    poly_s3_4 = zkey.poly_s3_4
    poly_qm = zkey.poly_qm
    poly_ql = zkey.poly_ql
    poly_qr = zkey.poly_qr
    poly_qo = zkey.poly_qo
    poly_qc = zkey.poly_qc
    poly_ls = zkey.poly_ls

    with log_elapsed_time("proof a, b, c", debug):
        tmp = Polynomial(Fr, coeffs=[b[1], b[0]]) * poly_zh_4
        tmp = tmp + poly_a
        poly_a_blinded_4 = Polynomial(Fr, coeffs=[b[1], b[0]]) * poly_zh_4 + poly_a
        poly_b_blinded_4 = Polynomial(Fr, coeffs=[b[3], b[2]]) * poly_zh_4 + poly_b
        poly_c_blinded_4 = Polynomial(Fr, coeffs=[b[5], b[4]]) * poly_zh_4 + poly_c

        proof_a = exp_tau(poly_a_blinded_4, taus)
        proof_b = exp_tau(poly_b_blinded_4, taus)
        proof_c = exp_tau(poly_c_blinded_4, taus)

    # round 2
    with log_elapsed_time("transcript1, 2", debug):
        transcript1 = b""
        transcript1 += int(bn128.normalize(zkey.qm)[0]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.qm)[1]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.ql)[0]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.ql)[1]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.qr)[0]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.qr)[1]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.qo)[0]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.qo)[1]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.qc)[0]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.qc)[1]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.s1)[0]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.s1)[1]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.s2)[0]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.s2)[1]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.s3)[0]).to_bytes(32, "big")
        transcript1 += int(bn128.normalize(zkey.s3)[1]).to_bytes(32, "big")
        for i in range(zkey.n_public):
            transcript1 += Fr.to_bytes(poly_a(omega**i))[::-1]
        transcript1 += int(proof_a[0]).to_bytes(32, "big")
        transcript1 += int(proof_a[1]).to_bytes(32, "big")
        transcript1 += int(proof_b[0]).to_bytes(32, "big")
        transcript1 += int(proof_b[1]).to_bytes(32, "big")
        transcript1 += int(proof_c[0]).to_bytes(32, "big")
        transcript1 += int(proof_c[1]).to_bytes(32, "big")
        h = keccak(transcript1)
        beta = Fr(int.from_bytes(h, "big"))
        transcript2 = Fr.to_bytes(beta)[::-1]
        h = keccak(transcript2)
        gamma = Fr(int.from_bytes(h, "big"))

    with log_elapsed_time("proof z", debug):
        with log_elapsed_time("prepare", debug, indent=1):
            evals = [Fr(1)]
            w = Fr(1)
            # poly_s1.calc_evals_if_necessary()
            # poly_s2.calc_evals_if_necessary()
            # poly_s3.calc_evals_if_necessary()
            omegas = [Fr(1)]
            for _ in range(zkey.domain_size):
                omegas.append(omegas[-1] * omega)
            assert omegas.pop() == Fr(1)
        with log_elapsed_time("poly", debug, indent=1):
            # tmp_an = poly_a + Polynomial(
            #     Fr, evals=[gamma + beta * omega for omega in omegas]
            # )
            # tmp_bn = poly_b + Polynomial(
            #     Fr, evals=[gamma + beta * zkey.k1 * omega for omega in omegas]
            # )
            # tmp_cn = poly_c + Polynomial(
            #     Fr, evals=[gamma + beta * zkey.k2 * omega for omega in omegas]
            # )
            # tmp_ad = (
            #     poly_a + Polynomial(Fr, evals=[gamma] * zkey.domain_size) + poly_s1 * beta
            # )
            # tmp_bd = (
            #     poly_b + Polynomial(Fr, evals=[gamma] * zkey.domain_size) + poly_s2 * beta
            # )
            # tmp_cd = (
            #     poly_c + Polynomial(Fr, evals=[gamma] * zkey.domain_size) + poly_s3 * beta
            # )
            # tmp_n = tmp_an * tmp_bn * tmp_cn
            # tmp_d = tmp_ad * tmp_bd * tmp_cd
            poly_a_evals = poly_a.evals
            poly_b_evals = poly_b.evals
            poly_c_evals = poly_c.evals
            poly_s1_evals = poly_s1.evals
            poly_s2_evals = poly_s2.evals
            poly_s3_evals = poly_s3.evals
            assert len(poly_a_evals) == zkey.domain_size
            assert len(poly_b_evals) == zkey.domain_size
            assert len(poly_c_evals) == zkey.domain_size
            assert len(poly_s1_evals) == zkey.domain_size
            assert len(poly_s2_evals) == zkey.domain_size
            assert len(poly_s3_evals) == zkey.domain_size
        with log_elapsed_time("evals", debug, indent=1):
            # tmp_n_evals = tmp_n.evals
            # tmp_d_evals = tmp_d.evals
            # assert len(tmp_n_evals) == zkey.domain_size
            # assert len(tmp_d_evals) == zkey.domain_size
            for i in range(zkey.domain_size):
                res = evals[-1]
                res = res * (poly_a_evals[i] + beta * w + gamma)
                res = res * (poly_b_evals[i] + beta * zkey.k1 * w + gamma)
                res = res * (poly_c_evals[i] + beta * zkey.k2 * w + gamma)
                res = res / (poly_a_evals[i] + poly_s1_evals[i] * beta + gamma)
                res = res / (poly_b_evals[i] + poly_s2_evals[i] * beta + gamma)
                res = res / (poly_c_evals[i] + poly_s3_evals[i] * beta + gamma)
                # res = res * tmp_n_evals[i] / tmp_d_evals[i]
                evals.append(res)
                w = w * omega
            assert evals.pop() == Fr(1)
        with log_elapsed_time("make poly", debug, indent=1):
            poly_z = Polynomial(Fr, evals=evals)
            poly_z_blinded_4 = (
                Polynomial(Fr, coeffs=[b[8], b[7], b[6]]) * poly_zh_4 + poly_z
            )
        with log_elapsed_time("proof z", debug, indent=1):
            proof_z = exp_tau(poly_z_blinded_4, taus)

    # round 3
    with log_elapsed_time("transcript3", debug):
        transcript3 = b""
        transcript3 += Fr.to_bytes(beta)[::-1]
        transcript3 += Fr.to_bytes(gamma)[::-1]
        transcript3 += int(proof_z[0]).to_bytes(32, "big")
        transcript3 += int(proof_z[1]).to_bytes(32, "big")
        h = keccak(transcript3)
        alpha = Fr(int.from_bytes(h, "big"))

    with log_elapsed_time("proof t", debug):
        with log_elapsed_time("pi", debug, indent=1):
            # pi_evals = [Fr(0) for _ in range(zkey.domain_size * 4)]
            # for i in range(zkey.domain_size * 4):
            #     for j in range(zkey.n_public):
            #         pi_evals[i] = pi_evals[i] - poly_l_4s[j].evals[i] * poly_a.evals[j]
            # poly_pi_4 = Polynomial(Fr, evals=pi_evals)
            poly_a_evals = poly_a.evals
            poly_pi_4 = Polynomial(Fr, evals=[0] * zkey.domain_size * 4)
            for j in range(zkey.n_public):
                poly_pi_4 = poly_pi_4 - poly_l_4s[j] * poly_a_evals[j]
        with log_elapsed_time("poly_z, poly_zw", debug, indent=1):
            poly_zw_blinded_4 = Polynomial(
                Fr,
                coeffs=[c * omega**i for i, c in enumerate(poly_z_blinded_4.coeffs)],
            )
            poly_z_blinded_8 = poly_z_blinded_4.extend(2 ** (zkey.power + 3))
            poly_zw_blinded_8 = poly_zw_blinded_4.extend(2 ** (zkey.power + 3))

        with log_elapsed_time("calc tzh", debug, indent=1):
            with log_elapsed_time("calc ta", debug, indent=2):
                ta = (
                    poly_a_blinded_4 * poly_b_blinded_4 * poly_qm_4
                    + poly_a_blinded_4 * poly_ql_4
                    + poly_b_blinded_4 * poly_qr_4
                    + poly_c_blinded_4 * poly_qo_4
                    + poly_pi_4
                    + poly_qc_4
                )
            with log_elapsed_time("calc tb", debug, indent=2):
                tb = (
                    (poly_a_blinded_4 + Polynomial(Fr, coeffs=[gamma, beta]))
                    * (poly_b_blinded_4 + Polynomial(Fr, coeffs=[gamma, beta * zkey.k1]))
                    * (poly_c_blinded_4 + Polynomial(Fr, coeffs=[gamma, beta * zkey.k2]))
                    * alpha
                    * poly_z_blinded_8
                )
            with log_elapsed_time("calc tc", debug, indent=2):
                tc = (
                    (poly_a_blinded_4 + poly_s1_4 * beta + gamma)
                    * (poly_b_blinded_4 + poly_s2_4 * beta + gamma)
                    * (poly_c_blinded_4 + poly_s3_4 * beta + gamma)
                    * alpha
                    * poly_zw_blinded_8
                )
            with log_elapsed_time("calc td", debug, indent=2):
                td = (poly_z_blinded_4 + Fr(-1)) * poly_l_4s[0] * alpha * alpha
            with log_elapsed_time("calc tzh", debug, indent=2):
                tzh = ta + td + (tb - tc)
                # tzh.calc_coeffs_if_necessary()

        with log_elapsed_time("t_coeffs", debug, indent=1):
            poly_t, tmp = divmod(tzh, SparsePolynomial(Fr, coeffs=[(0, -1), (zkey.domain_size, 1)]))
            t_coeffs = poly_t.coeffs
            assert all([int(c) == 0 for c in tmp.coeffs])
        # with log_elapsed_time("t_coeffs", debug, indent=1):
        #     t_coeffs = list(tzh.coeffs)
        #     for i in range(zkey.domain_size):
        #         t_coeffs[i] = t_coeffs[i] * Fr(-1)
        #     for i in range(zkey.domain_size, zkey.domain_size * 4 + 6):
        #         t_coeffs[i] = t_coeffs[i - zkey.domain_size] - t_coeffs[i]
        #         if i > zkey.domain_size * 3 + 5 and t_coeffs[i] != Fr(0):
        #             raise RuntimeError
        #     assert all(c == Fr(0) for c in t_coeffs[zkey.domain_size * 3 + 6 :])

        with log_elapsed_time("poly_t", debug, indent=1):
            poly_t_low = Polynomial(Fr, coeffs=t_coeffs[: zkey.domain_size] + [b[9]])
            poly_t_mid = Polynomial(
                Fr,
                coeffs=[t_coeffs[zkey.domain_size] - b[9]]
                + t_coeffs[zkey.domain_size + 1 : zkey.domain_size * 2]
                + [b[10]],
            )
            poly_t_high = Polynomial(
                Fr,
                coeffs=[t_coeffs[zkey.domain_size * 2] - b[10]]
                + t_coeffs[zkey.domain_size * 2 + 1 : zkey.domain_size * 3 + 6],
            )
        with log_elapsed_time("proof t", debug, indent=1):
            proof_t1 = exp_tau(poly_t_low, taus)
            proof_t2 = exp_tau(poly_t_mid, taus)
            proof_t3 = exp_tau(poly_t_high, taus)

    # round 4
    with log_elapsed_time("transcript4", debug):
        transcript4 = b""
        transcript4 += Fr.to_bytes(alpha)[::-1]
        transcript4 += int(proof_t1[0]).to_bytes(32, "big")
        transcript4 += int(proof_t1[1]).to_bytes(32, "big")
        transcript4 += int(proof_t2[0]).to_bytes(32, "big")
        transcript4 += int(proof_t2[1]).to_bytes(32, "big")
        transcript4 += int(proof_t3[0]).to_bytes(32, "big")
        transcript4 += int(proof_t3[1]).to_bytes(32, "big")
        h = keccak(transcript4)
        zeta = Fr(int.from_bytes(h, "big"))

    with log_elapsed_time("eval zeta", debug):
        eval_a = poly_a_blinded_4(zeta)
        eval_b = poly_b_blinded_4(zeta)
        eval_c = poly_c_blinded_4(zeta)
        eval_s1 = poly_s1(zeta)
        eval_s2 = poly_s2(zeta)
        eval_zw = poly_zw_blinded_4(zeta)

    # round 5
    with log_elapsed_time("transcript5", debug):
        transcript5 = b""
        transcript5 += Fr.to_bytes(zeta)[::-1]
        transcript5 += int(eval_a).to_bytes(32, "big")
        transcript5 += int(eval_b).to_bytes(32, "big")
        transcript5 += int(eval_c).to_bytes(32, "big")
        transcript5 += int(eval_s1).to_bytes(32, "big")
        transcript5 += int(eval_s2).to_bytes(32, "big")
        transcript5 += int(eval_zw).to_bytes(32, "big")
        h = keccak(transcript5)
        v = Fr(int.from_bytes(h, "big"))

    with log_elapsed_time("proof w_zeta, w_zetaw", debug):
        ra = (
            poly_qm * eval_a * eval_b
            + poly_ql * eval_a
            + poly_qr * eval_b
            + poly_qo * eval_c
            + poly_pi_4(zeta)
            + poly_qc
        )
        rb = (
            poly_z_blinded_4
            * (eval_a + beta * zeta + gamma)
            * (eval_b + beta * zkey.k1 * zeta + gamma)
            * (eval_c + beta * zkey.k2 * zeta + gamma)
            * alpha
        )
        rc = (
            (poly_s3 * beta + eval_c + gamma)
            * (eval_a + eval_s1 * beta + gamma)
            * (eval_b + eval_s2 * beta + gamma)
            * eval_zw
            * alpha
        )
        rd = (poly_z_blinded_4 + Fr(-1)) * poly_ls[0](zeta) * alpha * alpha
        re = (
            poly_t_low
            + poly_t_mid * zeta**zkey.domain_size
            + poly_t_high * zeta ** (zkey.domain_size * 2)
        ) * poly_zh_4(zeta)
        poly_r = ra + rb - rc + rd - re

        poly_w_zeta = (
            poly_r
            + (poly_a_blinded_4 - eval_a) * v
            + (poly_b_blinded_4 - eval_b) * v**2
            + (poly_c_blinded_4 - eval_c) * v**3
            + (poly_s1_4 - eval_s1) * v**4
            + (poly_s2_4 - eval_s2) * v**5
        )
        poly_w_zeta, tmp = divmod(poly_w_zeta, SparsePolynomial(Fr, coeffs=[(0, zeta * Fr(-1)), (1, 1)]))
        assert all([int(c) == 0 for c in tmp.coeffs])

        poly_w_zetaw = poly_z_blinded_4 - eval_zw
        poly_w_zetaw, tmp = divmod(poly_w_zetaw, SparsePolynomial(Fr, coeffs=[(0, zeta * omega * Fr(-1)), (1, 1)]))
        assert all([int(c) == 0 for c in tmp.coeffs])

        proof_w_zeta = exp_tau(poly_w_zeta, taus)
        proof_w_zetaw = exp_tau(poly_w_zetaw, taus)

    proof = Proof(
        a=proof_a + (FQ(1),),
        b=proof_b + (FQ(1),),
        c=proof_c + (FQ(1),),
        z=proof_z + (FQ(1),),
        t_low=proof_t1 + (FQ(1),),
        t_mid=proof_t2 + (FQ(1),),
        t_high=proof_t3 + (FQ(1),),
        w_zeta=proof_w_zeta + (FQ(1),),
        w_zetaw=proof_w_zetaw + (FQ(1),),
        eval_a=eval_a,
        eval_b=eval_b,
        eval_c=eval_c,
        eval_s1=eval_s1,
        eval_s2=eval_s2,
        eval_zw=eval_zw,
    )

    transcript6 = b""
    transcript6 += Fr.to_bytes(v)[::-1]
    transcript6 += int(proof_w_zeta[0]).to_bytes(32, "big")
    transcript6 += int(proof_w_zeta[1]).to_bytes(32, "big")
    transcript6 += int(proof_w_zetaw[0]).to_bytes(32, "big")
    transcript6 += int(proof_w_zetaw[1]).to_bytes(32, "big")
    h = keccak(transcript6)
    u = Fr(int.from_bytes(h, "big"))

    public_inputs = []
    for i in range(zkey.n_public):
        public_inputs.append(int(poly_a(omega**i)))

    if debug:
        print(f"{public_inputs = }")
        print(f"{beta = }")
        print(f"{gamma = }")
        print(f"{alpha = }")
        print(f"{zeta = }")
        print(f"{v = }")
        print(f"{u = }")

    return proof, public_inputs
