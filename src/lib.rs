mod hash;
mod field;
use crate::hash::keccak as _keccak;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

#[pyfunction]
fn keccak<'a>(py: Python<'a>, data: &'a [u8]) -> &'a PyBytes {
    let result = _keccak(data);
    PyBytes::new(py, &result)
}

fn register_hash_module(py: Python, parent_module: &PyModule) -> PyResult<()> {
    let hash_module = PyModule::new(py, "hash")?;
    hash_module.add_function(wrap_pyfunction!(keccak, hash_module)?)?;
    parent_module.add_submodule(hash_module)?;
    Ok(())
}

#[pymodule]
fn snarkpy(py: Python, m: &PyModule) -> PyResult<()> {
    register_hash_module(py, m)?;
    Ok(())
}
