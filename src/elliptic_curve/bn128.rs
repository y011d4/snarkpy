use num_bigint::BigInt;
use once_cell::sync::Lazy;
use pyo3::prelude::*;
use pyo3::types::PyType;

use crate::field::{
    integer_to_bigint, BigIntOrGFElement, GFElement, GFPolynomial, GFPolynomialElement, GF,
};
use crate::polynomial::SparsePolynomial;

pub static BN128: Lazy<BN128> = Lazy::new(|| BN128::new().unwrap());

#[pyclass]
#[derive(Clone)]
pub struct BN128 {
    #[pyo3(get)]
    pub fp: GF,
    #[pyo3(get)]
    pub fr: GF,
    #[pyo3(get)]
    pub b: GFElement,
    #[pyo3(get)]
    pub b2: GFPolynomialElement,
    #[pyo3(get)]
    pub b12: GFPolynomialElement,
    #[pyo3(get)]
    pub g: BN128Element,
    #[pyo3(get)]
    pub g2: BN128Element2,
    // pub g12: BN128Element12,
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct BN128Element {
    x: GFElement,
    y: GFElement,
    z: GFElement,
}

#[pyclass]
#[derive(Clone)]
pub struct BN128Element2 {
    x: GFPolynomialElement,
    y: GFPolynomialElement,
    z: GFPolynomialElement,
}

#[pyclass]
#[derive(Clone)]
pub struct BN128Element12 {
    x: GFPolynomialElement,
    y: GFPolynomialElement,
    z: GFPolynomialElement,
}

#[pymethods]
impl BN128 {
    #[new]
    pub fn new() -> PyResult<Self> {
        let fp = GF::new(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583"
                .parse()
                .unwrap(),
            "115792089237316195423570985008687907853269984665640564039457584007913129639936"
                .parse()
                .unwrap(),
        )?;
        let fr = GF::new(
            "21888242871839275222246405745257275088548364400416034343698204186575808495617"
                .parse()
                .unwrap(),
            "115792089237316195423570985008687907853269984665640564039457584007913129639936"
                .parse()
                .unwrap(),
        )?;
        let modulus2 = SparsePolynomial::new(
            fp.clone(),
            vec![
                (0, BigIntOrGFElement::BigInt("1".parse().unwrap())),
                (2, BigIntOrGFElement::BigInt("1".parse().unwrap())),
            ],
        );
        let modulus12 = SparsePolynomial::new(
            fp.clone(),
            vec![
                (0, BigIntOrGFElement::BigInt("2".parse().unwrap())),
                (6, BigIntOrGFElement::BigInt("-2".parse().unwrap())),
                (12, BigIntOrGFElement::BigInt("1".parse().unwrap())),
            ],
        );
        let fp2 = GFPolynomial::new(fp.clone(), modulus2)?;
        let fp12 = GFPolynomial::new(fp.clone(), modulus12)?;
        let b2 = fp2
            .__call__(vec![
                BigIntOrGFElement::BigInt("3".parse().unwrap()),
                BigIntOrGFElement::BigInt("0".parse().unwrap()),
            ])?
            .__truediv__(&fp2.__call__(vec![
                BigIntOrGFElement::BigInt("9".parse().unwrap()),
                BigIntOrGFElement::BigInt("1".parse().unwrap()),
            ])?)?;
        let b12 = fp12.__call__(vec![
            BigIntOrGFElement::BigInt("3".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
            BigIntOrGFElement::BigInt("0".parse().unwrap()),
        ])?;
        Ok(BN128 {
            fp: fp.clone(),
            fr,
            b: fp.__call__(BigIntOrGFElement::BigInt("3".parse().unwrap())),
            b2,
            b12,
            g: BN128Element {
                x: fp.__call__(BigIntOrGFElement::BigInt("1".parse().unwrap())),
                y: fp.__call__(BigIntOrGFElement::BigInt("2".parse().unwrap())),
                z: fp.__call__(BigIntOrGFElement::BigInt("1".parse().unwrap())),
            },
            g2: BN128Element2 {
                x: fp2.__call__(vec![
                        BigIntOrGFElement::BigInt("10857046999023057135944570762232829481370756359578518086990519993285655852781".parse().unwrap()),
                        BigIntOrGFElement::BigInt("11559732032986387107991004021392285783925812861821192530917403151452391805634".parse().unwrap()),
                ])?,
                y: fp2.__call__(vec![
                        BigIntOrGFElement::BigInt("8495653923123431417604973247489272438418190587263600148770280649306958101930".parse().unwrap()),
                        BigIntOrGFElement::BigInt("4082367875863433681332203403145435568316851327593401208105741076214120093531".parse().unwrap()),
                ])?,
                z: fp2.one(),
            },
            // g12:
        })
    }

    #[getter]
    fn p(&self) -> BigInt {
        integer_to_bigint(&BN128.fp.p)
    }

    #[getter]
    fn order(&self) -> BigInt {
        integer_to_bigint(&BN128.fr.p)
    }
}

#[pymethods]
impl BN128Element {
    #[new]
    pub fn new(x: BigIntOrGFElement, y: BigIntOrGFElement, z: Option<BigIntOrGFElement>) -> Self {
        Self {
            x: BN128.fp.__call__(x),
            y: BN128.fp.__call__(y),
            z: match z {
                Some(z) => BN128.fp.__call__(z),
                None => BN128
                    .fp
                    .__call__(BigIntOrGFElement::BigInt("1".parse().unwrap())),
            },
        }
    }

    fn __add__(&self, other: &Self) -> PyResult<Self> {
        if self.z == BN128.fp.zero() {
            return Ok(other.clone());
        }
        if other.z == BN128.fp.zero() {
            return Ok(self.clone());
        }
        let x1 = &self.x;
        let y1 = &self.y;
        let z1 = &self.z;
        let x2 = &other.x;
        let y2 = &other.y;
        let z2 = &other.z;
        let u1 = y2.__mul__(z1)?;
        let u2 = y1.__mul__(z2)?;
        let v1 = x2.__mul__(z1)?;
        let v2 = x1.__mul__(z2)?;
        if v1 == v2 && u1 == u2 {
            return Ok(self.double()?);
        } else if v1 == v2 {
            return Ok(Self {
                x: BN128.fp.one(),
                y: BN128.fp.one(),
                z: BN128.fp.zero(),
            });
        }
        let u = u1.__sub__(&u2)?;
        let v = v1.__sub__(&v2)?;
        let v_sq = v.__mul__(&v)?;
        let v_sq_v2 = v_sq.__mul__(&v2)?;
        let v_cu = v.__mul__(&v_sq)?;
        let w = z1.__mul__(&z2)?;
        let two = BN128
            .fp
            .__call__(BigIntOrGFElement::BigInt("2".parse().unwrap()));
        let a = u
            .__mul__(&u)?
            .__mul__(&w)?
            .__sub__(&v_cu)?
            .__sub__(&two.__mul__(&v_sq_v2)?)?;
        Ok(Self {
            x: v.__mul__(&a)?,
            y: u.__mul__(&v_sq_v2.__sub__(&a)?)?
                .__sub__(&v_cu.__mul__(&u2)?)?,
            z: v_cu.__mul__(&w)?,
        })
    }

    fn __mul__(&self, other: BigInt) -> PyResult<Self> {
        if &other == &BigInt::from(0) {
            return Ok(Self {
                x: BN128.fp.zero(),
                y: BN128.fp.one(),
                z: BN128.fp.zero(),
            });
        } else if &other == &BigInt::from(1) {
            return Ok(self.clone());
        } else if &other % 2 == BigInt::from(0) {
            return self.double()?.__mul__(&other / 2);
        } else {
            return self.double()?.__mul__(&other / 2)?.__add__(self);
        }
    }

    fn __eq__(&self, other: &Self) -> PyResult<bool> {
        let x1 = &self.x;
        let y1 = &self.y;
        let z1 = &self.z;
        let x2 = &other.x;
        let y2 = &other.y;
        let z2 = &other.z;
        let u1 = y2.__mul__(z1)?;
        let u2 = y1.__mul__(z2)?;
        let v1 = x2.__mul__(z1)?;
        let v2 = x1.__mul__(z2)?;
        Ok(u1 == u2 && v1 == v2)
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "BN128Element({}, {}, {})",
            self.x.__str__(),
            self.y.__str__(),
            self.z.__str__()
        ))
    }

    fn __str__(&self) -> PyResult<String> {
        let normalized = self.normalize()?;
        Ok(format!(
            "BN128Element({}, {})",
            normalized.x.__str__(),
            normalized.y.__str__(),
        ))
    }

    #[classmethod]
    fn zero(cls: &PyType) -> Self {
        Self {
            x: BN128.fp.zero(),
            y: BN128.fp.one(),
            z: BN128.fp.zero(),
        }
    }
}

impl BN128Element {
    fn double(&self) -> PyResult<Self> {
        let x = &self.x;
        let y = &self.y;
        let z = &self.z;
        let two = BN128
            .fp
            .__call__(BigIntOrGFElement::BigInt("2".parse().unwrap()));
        let three = BN128
            .fp
            .__call__(BigIntOrGFElement::BigInt("3".parse().unwrap()));
        let four = BN128
            .fp
            .__call__(BigIntOrGFElement::BigInt("4".parse().unwrap()));
        let eight = BN128
            .fp
            .__call__(BigIntOrGFElement::BigInt("8".parse().unwrap()));
        let w = three.__mul__(&x.__mul__(&x)?)?;
        let s = y.__mul__(&z)?;
        let b = x.__mul__(&y.__mul__(&s)?)?;
        let h = w.__mul__(&w)?.__sub__(&eight.__mul__(&b)?)?;
        let s_sq = s.__mul__(&s)?;
        Ok(Self {
            x: two.__mul__(&h)?.__mul__(&s)?,
            y: w.__mul__(&four.__mul__(&b)?.__sub__(&h)?)?
                .__sub__(&eight.__mul__(&y)?.__mul__(&y)?.__mul__(&s_sq)?)?,
            z: eight.__mul__(&s_sq.__mul__(&s)?)?,
        })
    }

    fn normalize(&self) -> PyResult<Self> {
        if self.z == BN128.fp.zero() {
            return Ok(self.clone());
        }
        let z_inv = self.z.__invert__();
        Ok(Self {
            x: self.x.__mul__(&z_inv)?,
            y: self.y.__mul__(&z_inv)?,
            z: BN128.fp.one(),
        })
    }
}
