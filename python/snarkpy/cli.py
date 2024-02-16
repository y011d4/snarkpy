import argparse
import json
from pathlib import Path

from snarkpy.plonk.prove import Proof, prove as plonk_prove
from snarkpy.plonk.verify import verify as plonk_verify


def prove(args):
    proof, public_inputs = plonk_prove(args.zkey, args.wtns, args.verbose)
    with args.proof.open("w") as fp:
        fp.write(proof.to_json())
    with args.public.open("w") as fp:
        json.dump([str(pi) for pi in public_inputs], fp)


def verify(args):
    with args.proof.open("r") as fp:
        proof = Proof.from_json(json.load(fp))
    with args.public.open("r") as fp:
        public = list(map(int, json.load(fp)))
    ok = plonk_verify(args.vkey, public, proof, args.verbose)
    if ok:
        print("Proof is valid")
        exit(0)
    else:
        print("Proof is invalid")
        exit(1)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    parser_plonk = subparsers.add_parser("plonk")
    subparsers_plonk = parser_plonk.add_subparsers()

    parser_plonk_prove = subparsers_plonk.add_parser("prove")
    parser_plonk_prove.add_argument("zkey", type=Path)
    parser_plonk_prove.add_argument("wtns", type=Path)
    parser_plonk_prove.add_argument("proof", type=Path)
    parser_plonk_prove.add_argument("public", type=Path)
    parser_plonk_prove.add_argument("-v", "--verbose", action="store_true")
    parser_plonk_prove.set_defaults(func=prove)

    parser_plonk_verify = subparsers_plonk.add_parser("verify")
    parser_plonk_verify.add_argument("vkey", type=Path)
    parser_plonk_verify.add_argument("public", type=Path)
    parser_plonk_verify.add_argument("proof", type=Path)
    parser_plonk_verify.add_argument("-v", "--verbose", action="store_true")
    parser_plonk_verify.set_defaults(func=verify)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
