# snarkpy

A zk-SNARK implemented in Python & Rust compatible with [snarkjs](https://github.com/iden3/snarkjs)

⚠️ It's not fully tested, so don't use it for production.

## Install

```sh
# Install maturin if not installed
pipx install maturin

# Install rustup if not installed
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

git clone https://github.com/y011d4/snarkpy.git
cd snarkpy
pytyhon3 -m venv venv
source venv/bin/activate
make install
```

## Usage

First, prepare `circuit_final.zkey`, `witness.wtns` and `verification_key.json` by referring to the [README in snarkjs](https://github.com/iden3/snarkjs).
If you only need test data, you can make them by just typing `make test-data`.

```sh
# plonk prove
snarkpy plonk prove zkey wtns proof.json public.json

# plonk verify
snarkpy plonk verify vkey public.json proof.json
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
