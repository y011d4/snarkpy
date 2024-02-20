use num_bigint::BigInt;
use num_traits::{FromBytes, One, Zero};
use once_cell::sync::Lazy;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};
use rug::{Assign, Integer};
use std::collections::HashMap;
use std::str::FromStr;
use std::sync::Mutex;

fn xgcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    let zero = BigInt::zero();
    if b == zero {
        (a, BigInt::one(), zero)
    } else {
        let (d, x, y) = xgcd(b.clone(), a.clone() % b.clone());
        (d, y.clone(), x - (a / b) * y)
    }
}

fn xgcd_integer(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    if b == &Integer::ZERO {
        (a.clone(), Integer::ONE.clone(), Integer::ZERO)
    } else {
        let c = Integer::from(a % b);
        let e = Integer::from(a / b);
        let (d, x, y) = xgcd_integer(&b, &c);
        (d, y.clone(), x - &e * &y)
    }
}

pub static NTH_ROOT_OF_UNITY_CACHE: Lazy<Mutex<HashMap<(Integer, Integer), GFElement>>> =
    Lazy::new(|| Mutex::new(HashMap::new()));

pub static R_INFO: Lazy<Mutex<HashMap<Integer, RInfo>>> = Lazy::new(|| Mutex::new(HashMap::new()));

#[derive(Clone)]
pub struct RInfo {
    r_bits: u64,
    r_mask: Integer,
    r2: Integer,
    r3: Integer,
    p_prime: Integer,
}

#[pyclass]
#[derive(Clone, Debug, PartialEq)]
pub struct GF {
    p: Integer,
    // r: BigInt,
    // r_bits: u64,
    // r_mask: Integer,
    // r2: Integer,
    // p_prime: Integer,
    // nth_root_of_unity_cache: HashMap<BigInt, GFElement>,
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct GFElement {
    #[pyo3(get)]
    gf: GF,
    value: Integer,
}

#[derive(Clone, Debug)]
pub enum BigIntOrGFElement {
    BigInt(BigInt),
    GFElement(GFElement),
}

#[pymethods]
impl GF {
    #[new]
    fn new(p: BigInt, r: BigInt) -> PyResult<Self> {
        // let integer_p = bigint_to_integer(&p);
        let p = bigint_to_integer(&p);
        let cache_key = &p;
        if R_INFO.lock().unwrap().contains_key(&cache_key) {
            return Ok(GF {
                p,
                // r_bits: info.r_bits,
                // r_mask: info.r_mask,
                // r2: info.r2,
                // p_prime: info.p_prime,
            });
        }
        let r_bits = r.bits() - 1;
        if &r & (&r - 1) != BigInt::from_str("0").unwrap() {
            return Err(PyValueError::new_err("r must be a power of 2"));
        }
        let r = bigint_to_integer(&r);
        if r <= p {
            return Err(PyValueError::new_err("r must be greater than p"));
        }
        let (mut g, mut p_inv) = (Integer::new(), Integer::new());
        (&mut g, &mut p_inv).assign(p.extended_gcd_ref(&r));
        if p_inv < Integer::ZERO {
            p_inv += &r;
        }
        // println!("g: {}, p_inv: {}", g, p_inv);
        // let (g, _, p_inv) = xgcd_integer(&r, &p);
        // if g != BigInt::from_str("1").unwrap() {
        //     return Err(PyValueError::new_err("p must be coprime with r"));
        // }
        if &g != Integer::ONE {
            return Err(PyValueError::new_err("p must be coprime with r"));
        }
        let r_1: Integer = Integer::from(&r - Integer::ONE);
        let r2: Integer = r.clone().pow_mod(&Integer::from(2).clone(), &p).unwrap();
        let r3: Integer = r.clone() * &r2 % &p;
        let info = RInfo {
            r_bits,
            r_mask: r_1,
            r2,
            r3,
            p_prime: -p_inv,
        };
        R_INFO
            .lock()
            .unwrap()
            .insert(cache_key.clone(), info.clone());
        Ok(GF {
            p,
            // r: r.clone(),
            // r_bits: r.bits() - 1,
            // r_mask: bigint_to_integer(&r_1),
            // r2: bigint_to_integer(&r.modpow(&BigInt::from_str("2").unwrap(), &p)),
            // p_prime: bigint_to_integer(&-p_inv),
            // nth_root_of_unity_cache: HashMap::new(),
        })
    }

    pub fn one(&self) -> GFElement {
        GFElement::new(BigInt::one(), self.clone(), false)
    }

    pub fn zero(&self) -> GFElement {
        GFElement::new(BigInt::zero(), self.clone(), false)
    }

    pub fn __call__(&self, value: BigIntOrGFElement) -> GFElement {
        match value {
            BigIntOrGFElement::BigInt(value) => GFElement::new(value, self.clone(), false),
            BigIntOrGFElement::GFElement(value) => value,
        }
    }

    fn __repr__(&self) -> String {
        format!("F_{}", self.p)
    }

    pub fn __str__(&self) -> String {
        format!("F_{}", self.p)
    }

    // pub fn __eq__(&self, other: Self) -> bool {
    //     self.p == other.p && self.r == other.r
    // }

    #[getter]
    fn n8(&self) -> u64 {
        // (integer_to_bits(&self.p).len() as u64 + 7) / 8
        (self.p.signed_bits() - 1) as u64
    }

    fn from_bytes(&self, data: &[u8], is_montgomery: bool) -> GFElement {
        let value = BigInt::from_bytes_le(num_bigint::Sign::Plus, data);
        // value %= &self.p.try_into().unwrap();
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
        value.gf == self.clone()
    }

    pub fn nth_root_of_unity(&self, n: BigInt) -> PyResult<GFElement> {
        if n.sign() == num_bigint::Sign::Minus {
            return Err(PyValueError::new_err("n must be positive"));
        }
        // if self.nth_root_of_unity_cache.contains_key(&n) {
        let n: Integer = bigint_to_integer(&n);
        let cache_key = (self.p.clone(), n.clone());
        if NTH_ROOT_OF_UNITY_CACHE
            .lock()
            .unwrap()
            .contains_key(&cache_key)
        {
            // return Ok(self.nth_root_of_unity_cache[&n].clone());
            return Ok(NTH_ROOT_OF_UNITY_CACHE.lock().unwrap()[&cache_key].clone());
        }
        let p_1 = Integer::from(&self.p - Integer::ONE);
        let tmp = Integer::from(&p_1 % &n);
        if !tmp.is_zero() {
            return Err(PyValueError::new_err("n must divide p - 1"));
        }
        let k = self.__call__(BigIntOrGFElement::BigInt(BigInt::from_str("5").unwrap()));
        let res = GFElement {
            value: self.pow(&k.value, &(&p_1 / n)),
            gf: self.clone(),
        };
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
            GFElement {
                gf,
                value: bigint_to_integer(&value),
            }
        } else {
            GFElement {
                gf: gf.clone(),
                value: gf.to_montgomery(&bigint_to_integer(&value)),
            }
        }
    }

    pub fn __add__(&self, other: &Self) -> PyResult<Self> {
        // if self.gf != other.gf {
        //     return Err(PyValueError::new_err("GF mismatch"));
        // }
        Ok(GFElement {
            gf: other.gf.clone(),
            value: self.gf.add(&self.value, &other.value),
        })
    }

    pub fn __sub__(&self, other: &Self) -> PyResult<Self> {
        // if self.gf != other.gf {
        //     return Err(PyValueError::new_err("GF mismatch"));
        // }
        Ok(GFElement {
            value: self.gf.sub(&self.value, &other.value),
            gf: other.gf.clone(),
        })
    }

    pub fn __mul__(&self, other: &Self) -> PyResult<Self> {
        // if self.gf != other.gf {
        //     return Err(PyValueError::new_err("GF mismatch"));
        // }
        let tmp = Integer::from(&other.value * &self.value);
        Ok(GFElement {
            value: self.gf.from_montgomery(&tmp),
            gf: other.gf.clone(),
        })
    }

    fn __invert__(&self) -> Self {
        GFElement {
            value: self.gf.inv(&self.value),
            gf: self.gf.clone(),
        }
    }

    pub fn __truediv__(&self, other: &Self) -> PyResult<Self> {
        // if self.gf != other.gf {
        //     return Err(PyValueError::new_err("GF mismatch"));
        // }
        Ok(GFElement {
            value: self.gf.div(&self.value, &other.value),
            gf: other.gf.clone(),
        })
    }

    pub fn __pow__(&self, other: BigInt, p: Option<BigInt>) -> PyResult<Self> {
        if let Some(p) = p {
            let p: Integer = bigint_to_integer(&p);
            if p != self.gf.p {
                return Err(PyValueError::new_err("GF mismatch"));
            }
        }
        Ok(GFElement {
            value: self.gf.pow(&self.value, &bigint_to_integer(&other)),
            gf: self.gf.clone(),
        })
    }

    fn __int__(&self) -> BigInt {
        integer_to_bigint(&self.gf.from_montgomery(&self.value))
    }

    pub fn __eq__(&self, other: Self) -> bool {
        self.gf == other.gf && self.value == other.value
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

impl<'a> FromPyObject<'a> for BigIntOrGFElement {
    fn extract(ob: &'a PyAny) -> PyResult<Self> {
        if let Ok(s) = BigInt::extract(ob) {
            Ok(BigIntOrGFElement::BigInt(s))
        } else if let Ok(s) = GFElement::extract(ob) {
            Ok(BigIntOrGFElement::GFElement(s))
        } else {
            Err(PyValueError::new_err("Invalid type"))
        }
    }
}

impl GF {
    fn add(&self, a: &Integer, b: &Integer) -> Integer {
        let tmp = Integer::from(a + b);
        if tmp >= self.p {
            tmp - &self.p
        } else {
            tmp
        }
    }

    fn sub(&self, a: &Integer, b: &Integer) -> Integer {
        let c = a.clone() - b;
        let result = c % &self.p;
        if &result < &Integer::ZERO {
            result + &self.p
        } else {
            result
        }
    }

    fn mul(&self, a: &Integer, b: &Integer) -> Integer {
        let tmp = Integer::from(a * b);
        self.from_montgomery(&tmp)
    }

    fn inv(&self, a: &Integer) -> Integer {
        let mut a_inv = a.clone();
        a_inv.invert_mut(&self.p).unwrap();
        // let mut a_inv = xgcd_integer(a, &self.p).1;
        if a_inv < 0 {
            a_inv += &self.p;
        }
        let r3 = R_INFO.lock().unwrap()[&self.p].r3.clone();
        self.from_montgomery(&(a_inv * &r3))
        // self.pow(a, &(&self.p - Integer::from(2)))
    }

    fn div(&self, a: &Integer, b: &Integer) -> Integer {
        self.mul(a, &self.inv(b))
    }

    fn pow(&self, a: &Integer, b: &Integer) -> Integer {
        let p_1 = Integer::from(&self.p - Integer::ONE);
        let mut b = Integer::from(b % &p_1);
        // if &b.clone().signum() == Integer::NEG_ONE {
        //     b += &p_1;
        // }
        if b < Integer::ZERO {
            b += &p_1;
        }
        // let b_string = b.to_str_radix(2);
        // let bits = b_string.chars();
        //
        let bitnum = b.signed_bits();
        let mut a0 = self.to_montgomery(Integer::ONE);
        let mut a1 = a.clone();
        for i in (0..bitnum - 1).rev() {
            let bit_is_one = b.get_bit(i);
            if bit_is_one {
                a0 = self.mul(&a0, &a1);
                a1 = self.mul(&a1, &a1);
            } else {
                a1 = self.mul(&a0, &a1);
                a0 = self.mul(&a0, &a0);
            }
        }
        //
        // let b_bits = integer_to_bits(&b);
        // let mut a0 = self.to_montgomery(Integer::ONE);
        // let mut a1 = a.clone();
        // for bit in b_bits.iter().rev() {
        //     match bit {
        //         0 => {
        //             a1 = self.mul(&a0, &a1);
        //             a0 = self.mul(&a0, &a0);
        //         }
        //         1 => {
        //             a0 = self.mul(&a0, &a1);
        //             a1 = self.mul(&a1, &a1);
        //         }
        //         _ => panic!("Invalid bit"),
        //     }
        // }
        a0
    }

    fn to_montgomery(&self, value: &Integer) -> Integer {
        let tmp = value.clone() * &R_INFO.lock().unwrap()[&self.p].r2;
        self.from_montgomery(&tmp)
    }

    fn from_montgomery(&self, value: &Integer) -> Integer {
        let r_info = &R_INFO.lock().unwrap()[&self.p];
        let mut tmp = Integer::from(value * &r_info.p_prime);
        tmp &= &r_info.r_mask;
        tmp *= &self.p;
        tmp += value;
        tmp >>= r_info.r_bits as usize;
        tmp.shrink_to_fit();
        if tmp >= self.p {
            tmp - &self.p
        } else {
            tmp
        }
    }
}

fn integer_to_bigint(value: &Integer) -> BigInt {
    let mut value = value.clone();
    let mut sign = 1;
    if value < Integer::ZERO {
        value = -value;
        sign = -1;
    }
    // let mut tmp = Integer::ONE.clone();
    // let tmp_256: Integer = "256".parse().unwrap();
    // let mut n = 1;
    // while tmp <= value {
    //     tmp *= &tmp_256;
    //     n += 1;
    // }
    // let mut digits = vec![0; n];
    // value.write_digits(&mut digits, Order::Lsf);
    let mut digits = integer_to_bytes(&value);
    digits.push(0); // BigInt::from_bytes_le が負にならないように
    let mut result = BigInt::from_le_bytes(&digits);
    if sign == -1 {
        result = -result;
    }
    result
}

fn integer_to_bytes(value: &Integer) -> Vec<u8> {
    let u8_mask = Integer::from(0xff);
    let mut result = Vec::new();
    let mut value = value.clone();
    while value > Integer::ZERO {
        let tmp = Integer::from(&value & &u8_mask);
        result.push(tmp.to_u8().unwrap());
        value >>= 8;
    }
    result
}

fn bigint_to_integer(value: &BigInt) -> Integer {
    value.to_str_radix(10).parse().unwrap()
}

fn integer_to_bits(value: &Integer) -> Vec<u8> {
    let mut result = Vec::new();
    let mut value = value.clone();
    while value > Integer::ZERO {
        // let tmp = Integer::from(&value & Integer::ONE);
        // result.push(tmp.to_u8().unwrap());
        let tmp = value.get_bit(0);
        result.push(tmp as u8);
        value >>= 1;
    }
    result
}
