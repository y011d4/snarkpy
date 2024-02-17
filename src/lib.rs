mod field;
mod hash;
mod polynomial;
use crate::field::{GFElement, GF};
use crate::hash::keccak as _keccak;
use crate::polynomial::Polynomial;
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

fn register_field_module(py: Python, parent_module: &PyModule) -> PyResult<()> {
    let field_module = PyModule::new(py, "field")?;
    field_module.add_class::<GF>()?;
    field_module.add_class::<GFElement>()?;
    parent_module.add_submodule(field_module)?;
    Ok(())
}

fn register_polynomial_module(py: Python, parent_module: &PyModule) -> PyResult<()> {
    let polynomial_module = PyModule::new(py, "polynomial")?;
    polynomial_module.add_class::<Polynomial>()?;
    parent_module.add_submodule(polynomial_module)?;
    Ok(())
}

#[pymodule]
fn _snarkpy(py: Python, m: &PyModule) -> PyResult<()> {
    register_hash_module(py, m)?;
    register_field_module(py, m)?;
    register_polynomial_module(py, m)?;
    Ok(())
}
