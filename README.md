# snarkpy

A zk-SNARK implemented in Python & Rust compatible with [snarkjs](https://github.com/iden3/snarkjs)

⚠️ It's not fully tested, so don't use it for production.

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

## Usage

First, prepare `circuit_final.zkey`, `witness.wtns` and `verification_key.json` by referring to the [README in snarkjs](https://github.com/iden3/snarkjs).

```sh
# plonk prove
snarkpy plonk prove circuit_final.zkey witness.wtns proof.json public.json

# plonk verify
snarkpy plonk verify circuit_final.zkey witness.wtns proof.json public.json
```

## TODO

- Implement other algorithms than PLONK.
- Implement EC by myself and replace py_ecc.
- Replace mathematical implementation like Polynomial or GF with implementation Rust using PyO3.

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
