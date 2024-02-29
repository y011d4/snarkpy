FROM python:3.11-buster as builder
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs > /tmp/install_rust.sh && sh /tmp/install_rust.sh -y && rm /tmp/install_rust.sh
RUN pip3 install maturin
WORKDIR /opt/app
COPY pyproject.toml /opt/app/pyproject.toml
COPY Makefile /opt/app/Makefile
COPY README.md /opt/app/README.md
COPY Cargo.toml /opt/app/Cargo.toml
COPY python /opt/app/python
COPY src /opt/app/src
# RUN bash -c "source $HOME/.cargo/env && pip install -r requirements.txt && python setup.py install"
RUN bash -c "source $HOME/.cargo/env && make install"

FROM python:3.11-slim-buster as runner
COPY --from=builder /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=builder /usr/local/bin/snarkpy /usr/local/bin/snarkpy
CMD python
