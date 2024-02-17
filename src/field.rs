use num_bigint::BigInt;
use once_cell::sync::Lazy;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};
use std::collections::HashMap;
use std::str::FromStr;
use std::sync::Mutex;

fn xgcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    let zero = BigInt::from_str("0").unwrap();
    if b == zero {
        (a, BigInt::from_str("1").unwrap(), zero)
    } else {
        let (d, x, y) = xgcd(b.clone(), a.clone() % b.clone());
        (d, y.clone(), x - (a / b) * y)
    }
}

pub static NTH_ROOT_OF_UNITY_CACHE: Lazy<Mutex<HashMap<(BigInt, BigInt), GFElement>>> =
    Lazy::new(|| Mutex::new(HashMap::new()));

#[pyclass]
#[derive(Clone, Debug)]
pub struct GF {
    #[pyo3(get)]
    p: BigInt,
    r: BigInt,
    r_bits: u64,
    r_mask: BigInt,
    r2: BigInt,
    p_prime: BigInt,
    // nth_root_of_unity_cache: HashMap<BigInt, GFElement>,
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct GFElement {
    #[pyo3(get)]
    gf: GF,
    value: BigInt,
}

#[derive(Clone, Debug)]
pub enum Value {
    BigInt(BigInt),
    GFElement(GFElement),
}

#[pymethods]
impl GF {
    #[new]
    fn new(p: BigInt, r: BigInt) -> PyResult<Self> {
        if &r & (&r - 1) != BigInt::from_str("0").unwrap() {
            return Err(PyValueError::new_err("r must be a power of 2"));
        }
        if r <= p {
            return Err(PyValueError::new_err("r must be greater than p"));
        }
        let (g, _, p_inv) = xgcd(r.clone(), p.clone());
        if g != BigInt::from_str("1").unwrap() {
            return Err(PyValueError::new_err("p must be coprime with r"));
        }
        Ok(GF {
            p: p.clone(),
            r: r.clone(),
            r_bits: r.bits() - 1,
            r_mask: &r - 1,
            r2: r.modpow(&BigInt::from_str("2").unwrap(), &p),
            p_prime: -p_inv,
            // nth_root_of_unity_cache: HashMap::new(),
        })
    }

    fn add(&self, a: BigInt, b: BigInt) -> BigInt {
        (a + b) % &self.p
    }

    fn sub(&self, a: BigInt, b: BigInt) -> BigInt {
        let result = (a - b) % &self.p;
        if result.sign() == num_bigint::Sign::Minus {
            result + &self.p
        } else {
            result
        }
    }

    fn mul(&self, a: BigInt, b: BigInt) -> BigInt {
        self.from_montgomery(a * b)
    }

    fn inv(&self, a: BigInt) -> BigInt {
        self.pow(a, &self.p - BigInt::from_str("2").unwrap())
    }

    fn div(&self, a: BigInt, b: BigInt) -> BigInt {
        self.mul(a, self.inv(b))
    }

    fn pow(&self, a: BigInt, b: BigInt) -> BigInt {
        let p_1 = &self.p - BigInt::from_str("1").unwrap();
        let mut b = b % &p_1;
        if b.sign() == num_bigint::Sign::Minus {
            b += &p_1;
        }
        let b_string = b.to_str_radix(2);
        let bits = b_string.chars();
        let mut a0 = self.to_montgomery(BigInt::from_str("1").unwrap());
        let mut a1 = a;
        for bit in bits {
            match bit {
                '0' => {
                    a1 = self.mul(a0.clone(), a1.clone());
                    a0 = self.mul(a0.clone(), a0.clone());
                }
                '1' => {
                    a0 = self.mul(a0.clone(), a1.clone());
                    a1 = self.mul(a1.clone(), a1.clone());
                }
                _ => panic!("Invalid bit"),
            }
        }
        a0
    }

    pub fn one(&self) -> GFElement {
        GFElement::new(BigInt::from_str("1").unwrap(), self.clone(), false)
    }

    pub fn zero(&self) -> GFElement {
        GFElement::new(BigInt::from_str("0").unwrap(), self.clone(), false)
    }

    pub fn __call__(&self, value: Value) -> GFElement {
        match value {
            Value::BigInt(value) => GFElement::new(value, self.clone(), false),
            Value::GFElement(value) => value,
        }
    }

    fn to_montgomery(&self, value: BigInt) -> BigInt {
        self.from_montgomery(value * &self.r2)
    }

    fn from_montgomery(&self, value: BigInt) -> BigInt {
        let t = (value.clone() + ((value.clone() * &self.p_prime) & &self.r_mask) * &self.p)
            >> self.r_bits;
        if t >= self.p {
            t - &self.p
        } else {
            t - 0
        }
    }

    fn __repr__(&self) -> String {
        format!("F_{}", self.p)
    }

    pub fn __str__(&self) -> String {
        format!("F_{}", self.p)
    }

    pub fn __eq__(&self, other: Self) -> bool {
        self.p == other.p && self.r == other.r
    }

    #[getter]
    fn n8(&self) -> u64 {
        (self.p.bits() + 7) / 8
    }

    fn from_bytes(&self, data: &[u8], is_montgomery: bool) -> GFElement {
        let mut value = BigInt::from_bytes_le(num_bigint::Sign::Plus, &data);
        value = value % &self.p;
        GFElement::new(value, self.clone(), is_montgomery)
    }

    fn to_bytes<'a>(&self, value: GFElement, py: Python<'a>) -> &'a PyBytes {
        let value = value.__int__();
        let (_, mut result) = value.to_bytes_le();
        let n8 = self.n8();
        if result.len() < n8 as usize {
            result.resize(n8 as usize, 0);
        }
        PyBytes::new(py, &result)
    }

    fn __contains__(&self, value: GFElement) -> bool {
        value.gf.__eq__(self.clone())
    }

    pub fn nth_root_of_unity(&mut self, n: BigInt) -> PyResult<GFElement> {
        if n.sign() == num_bigint::Sign::Minus {
            return Err(PyValueError::new_err("n must be positive"));
        }
        // if self.nth_root_of_unity_cache.contains_key(&n) {
        let cache_key = (self.p.clone(), n.clone());
        if NTH_ROOT_OF_UNITY_CACHE
            .lock()
            .unwrap()
            .contains_key(&cache_key)
        {
            // return Ok(self.nth_root_of_unity_cache[&n].clone());
            return Ok(NTH_ROOT_OF_UNITY_CACHE.lock().unwrap()[&cache_key].clone());
        }
        let p_1 = &self.p - BigInt::from_str("1").unwrap();
        if &p_1 % &n != BigInt::from_str("0").unwrap() {
            return Err(PyValueError::new_err("n must divide p - 1"));
        }
        let k = self.__call__(Value::BigInt(BigInt::from_str("5").unwrap()));
        let res = k.__pow__(&p_1 / n.clone(), None)?;
        NTH_ROOT_OF_UNITY_CACHE
            .lock()
            .unwrap()
            .insert(cache_key, res.clone());
        Ok(res)
    }
}

#[pymethods]
impl GFElement {
    #[new]
    pub fn new(value: BigInt, gf: GF, is_montgomery: bool) -> Self {
        if is_montgomery {
            GFElement { gf, value }
        } else {
            GFElement {
                gf: gf.clone(),
                value: gf.to_montgomery(value),
            }
        }
    }

    pub fn __add__(&self, other: Self) -> PyResult<Self> {
        if !self.gf.__eq__(other.gf) {
            return Err(PyValueError::new_err("GF mismatch"));
        }
        Ok(GFElement::new(
            self.gf.add(self.value.clone(), other.value),
            self.gf.clone(),
            true,
        ))
    }

    pub fn __sub__(&self, other: Self) -> PyResult<Self> {
        if !self.gf.__eq__(other.gf) {
            return Err(PyValueError::new_err("GF mismatch"));
        }
        Ok(GFElement::new(
            self.gf.sub(self.value.clone(), other.value),
            self.gf.clone(),
            true,
        ))
    }

    pub fn __mul__(&self, other: Self) -> PyResult<Self> {
        if !self.gf.__eq__(other.gf) {
            return Err(PyValueError::new_err("GF mismatch"));
        }
        Ok(GFElement::new(
            self.gf.mul(self.value.clone(), other.value),
            self.gf.clone(),
            true,
        ))
    }

    fn __invert__(&self) -> Self {
        GFElement::new(self.gf.inv(self.value.clone()), self.gf.clone(), true)
    }

    pub fn __truediv__(&self, other: Self) -> PyResult<Self> {
        if !self.gf.__eq__(other.gf) {
            return Err(PyValueError::new_err("GF mismatch"));
        }
        Ok(GFElement::new(
            self.gf.div(self.value.clone(), other.value),
            self.gf.clone(),
            true,
        ))
    }

    pub fn __pow__(&self, other: BigInt, p: Option<BigInt>) -> PyResult<Self> {
        match p {
            Some(p) => {
                if p != self.gf.p {
                    return Err(PyValueError::new_err("GF mismatch"));
                }
            }
            None => {}
        }
        Ok(GFElement::new(
            self.gf.pow(self.value.clone(), other),
            self.gf.clone(),
            true,
        ))
    }

    fn __int__(&self) -> BigInt {
        self.gf.from_montgomery(self.value.clone())
    }

    pub fn __eq__(&self, other: Self) -> bool {
        self.gf.__eq__(other.gf) && self.value == other.value
    }

    fn __repr__(&self) -> String {
        format!("{} in {}", self.__int__(), self.gf.__str__())
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.__int__())
    }

    fn __deepcopy__(&self, _memo: &PyDict) -> Self {
        self.clone()
    }
}

impl<'a> FromPyObject<'a> for Value {
    fn extract(ob: &'a PyAny) -> PyResult<Self> {
        if let Ok(s) = BigInt::extract(ob) {
            Ok(Value::BigInt(s))
        } else if let Ok(s) = GFElement::extract(ob) {
            Ok(Value::GFElement(s))
        } else {
            Err(PyValueError::new_err("Invalid type"))
        }
    }
}
