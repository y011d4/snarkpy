# snarkpy

A zk-SNARK implemented in Python compatible with [snarkjs](https://github.com/iden3/snarkjs)

## Install

```sh
# Install poetry if not installed
curl -sSL https://install.python-poetry.org | python3 -

# Install rustup if not installed
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

git clone https://github.com/y011d4/snarkpy.git
poetry install
poetry run maturin develop
```

### Usage

First, prepare `circuit_final.zkey`, `witness.wtns` and `verification_key.json` by referring to the [README in snarkjs](https://github.com/iden3/snarkjs).

```sh
# plonk prove
snarkpy plonk prove circuit_final.zkey witness.wtns proof.json public.json

# plonk verify
snarkpy plonk verify circuit_final.zkey witness.wtns proof.json public.json
```

### TODO

- Implement other algorithms than PLONK.
- Implement EC by myself and replace py_ecc.
- Replace mathematical implementation like Polynomial or GF with implementation Rust using PyO3.
