.PHONY: all
all: install

.PHONY: install
install:
	maturin develop --release

.PHONY: install-dev
install-dev:
	maturin develop --release -E dev

.PHONY: test
test: install-dev
	pytest tests -vv

.PHONY: test-data
test-data: zkey vkey wtns

circom:
	wget https://github.com/iden3/circom/releases/latest/download/circom-linux-amd64 -O circom
	chmod u+x circom

ptau:
	wget https://storage.googleapis.com/zkevm/ptau/powersOfTau28_hez_final_12.ptau -O ptau

circuits: circom
	npm install
	mkdir circuits
	./circom circuit.circom --wasm --r1cs --sym -o circuits

zkey: circuits ptau
	npx snarkjs plonk setup circuits/circuit.r1cs ptau zkey

vkey: circuits zkey
	npx snarkjs zkey export verificationkey zkey vkey

wtns: circuits
	echo '{ "inputs": 1337 }' > circuits/input.json
	node circuits/circuit_js/generate_witness.js circuits/circuit_js/circuit.wasm circuits/input.json wtns
