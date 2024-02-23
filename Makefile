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
	pytest tests
