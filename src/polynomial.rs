use num_bigint::BigInt;
use once_cell::sync::Lazy;
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};
use rayon::prelude::*;
use std::str::FromStr;

use crate::field::{BigIntOrGFElement, GFElement, GF, OMEGA_POWS};

static MAX_K: u32 = 20; // TODO: 28
static MAX_L: u32 = 2u32.pow(MAX_K);
static REVBITS: Lazy<Vec<u32>> = Lazy::new(|| {
    fn reverse_bits(n: u32) -> u32 {
        let mut n = (n >> 16) | (n << 16);
        n = ((n & 0xff00ff00) >> 8) | ((n & 0x00ff00ff) << 8);
        n = ((n & 0xf0f0f0f0) >> 4) | ((n & 0x0f0f0f0f) << 4);
        n = ((n & 0xcccccccc) >> 2) | ((n & 0x33333333) << 2);
        n = ((n & 0xaaaaaaaa) >> 1) | ((n & 0x55555555) << 1);
        n
    }
    let mut revbits = Vec::new();
    for i in 0..MAX_L {
        revbits.push(reverse_bits(i) >> (32 - MAX_K));
    }
    revbits
});

#[pyclass]
#[derive(Clone, Debug)]
pub struct Polynomial {
    coeffs: Option<Vec<GFElement>>,
    evals: Option<Vec<GFElement>>,
    #[pyo3(get)]
    gf: GF,
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct SparsePolynomial {
    #[pyo3(get)]
    coeffs: Vec<(usize, GFElement)>,
    // evals: Option<Vec<GFElement>>,
    #[pyo3(get)]
    gf: GF,
}

#[derive(Clone, Debug)]
pub enum PolynomialOrGFElement {
    Polynomial(Polynomial),
    GFElement(GFElement),
}

#[pymethods]
impl Polynomial {
    #[new]
    fn new(
        gf: GF,
        coeffs: Option<Vec<BigIntOrGFElement>>,
        evals: Option<Vec<BigIntOrGFElement>>,
    ) -> PyResult<Self> {
        if coeffs.is_none() && evals.is_none() {
            return Err(PyValueError::new_err(
                "coeffs and evals cannot both be None",
            ));
        }
        let coeffs = match coeffs {
            Some(coeffs) => {
                let l = coeffs.len();
                let n = match l {
                    1 => 1,
                    _ => 2usize.pow((l - 1).ilog2() + 1),
                };
                let mut result = Vec::new();
                for coeff in coeffs {
                    match coeff {
                        BigIntOrGFElement::BigInt(value) => {
                            result.push(GFElement::new(value, gf.clone(), false));
                        }
                        BigIntOrGFElement::GFElement(value) => {
                            result.push(value);
                        }
                    }
                }
                for _ in 0..(n - l) {
                    result.push(gf.zero());
                }
                Some(result)
            }
            None => None,
        };
        let evals = match evals {
            Some(evals) => {
                let l = evals.len();
                let n = 2usize.pow((l - 1).ilog2() + 1);
                let mut result = Vec::new();
                for eval in evals {
                    match eval {
                        BigIntOrGFElement::BigInt(value) => {
                            result.push(GFElement::new(value, gf.clone(), false));
                        }
                        BigIntOrGFElement::GFElement(value) => {
                            result.push(value);
                        }
                    }
                }
                for _ in 0..(n - l) {
                    result.push(gf.zero());
                }
                Some(result)
            }
            None => None,
        };
        Ok(Polynomial { coeffs, evals, gf })
    }

    fn __call__(&mut self, x: BigIntOrGFElement) -> PyResult<GFElement> {
        self.calc_coeffs_if_necessary()?;
        let x = match x {
            BigIntOrGFElement::BigInt(value) => GFElement::new(value, self.gf.clone(), false),
            BigIntOrGFElement::GFElement(value) => value,
        };
        let mut result = self.gf.zero();
        for c in self.coeffs.as_ref().unwrap().iter().rev() {
            result = result.__mul__(&x).unwrap().__add__(&c).unwrap();
        }
        Ok(result)
    }

    // fn fft(&self, coeffs: Vec<GFElement>, omegas: Vec<GFElement>) -> PyResult<Vec<GFElement>> {
    //     let n = coeffs.len();
    //     if n & (n - 1) != 0 {
    //         return Err(PyValueError::new_err("n must be a power of 2"));
    //     }
    //     if n == 1 {
    //         return Ok(coeffs);
    //     }
    //     let ye = self.fft(coeffs.iter().step_by(2).cloned().collect(), omegas.clone())?;
    //     let yo = self.fft(
    //         coeffs.iter().skip(1).step_by(2).cloned().collect(),
    //         omegas.clone(),
    //     )?;
    //     let omega = &omegas[n.ilog2() as usize];
    //     let mut omega_i = self.gf.__call__(Value::GFElement(self.gf.one()));
    //     let mut evals = vec![self.gf.zero(); n];
    //     for i in 0..(n / 2) {
    //         let tmp = omega_i.__mul__(yo[i].clone())?;
    //         evals[i] = ye[i].__add__(tmp.clone())?;
    //         evals[i + n / 2] = ye[i].__sub__(tmp)?;
    //         omega_i = omega_i.__mul__(omega.clone())?;
    //     }
    //     Ok(evals)
    // }

    // fn fft(coeffs: Vec<GFElement>, omegas: Vec<GFElement>) -> Vec<GFElement> {
    //     let n = coeffs.len();
    //     if n == 1 {
    //         return coeffs;
    //     }
    //     let mut even = Vec::new();
    //     let mut odd = Vec::new();
    //     for i in 0..n {
    //         if i % 2 == 0 {
    //             even.push(coeffs[i]);
    //         } else {
    //             odd.push(coeffs[i]);
    //         }
    //     }
    //     let mut even_fft = fft(even, omegas.iter().step_by(2).cloned().collect());
    //     let mut odd_fft = fft(odd, omegas.iter().step_by(2).cloned().collect());
    //     let mut result =
    //         vec![GFElement::new(BigInt::from_str("0").unwrap(), self.gf.clone(), false); n];
    //     for i in 0..n / 2 {
    //         let t = omegas[i] * odd_fft[i];
    //         result[i] = even_fft[i] + t;
    //         result[i + n / 2] = even_fft[i] - t;
    //     }
    //     result
    // }

    fn __repr__(&mut self) -> PyResult<String> {
        self.calc_coeffs_if_necessary()?;
        let mut res = Vec::new();
        for (i, c) in self.coeffs.as_ref().unwrap().iter().enumerate() {
            if c.__eq__(self.gf.zero()) {
                continue;
            }
            if i == 0 {
                res.push(c.__str__());
            } else if c.__eq__(self.gf.one()) {
                res.push(format!("x^{}", i));
            } else {
                res.push(format!("{} * x^{}", c.__str__(), i));
            }
        }
        let mut ret = res.join(" + ");
        if ret.len() > 1024 {
            ret = ret[..1024].to_string() + " ...";
        }
        Ok(format!("{} in {}", ret, self.gf.__str__()))
    }

    fn __len__(&self) -> usize {
        let mut len_coeffs = 0;
        let mut len_evals = 0;
        if let Some(ref coeffs) = self.coeffs {
            len_coeffs = coeffs.len();
        }
        if let Some(ref evals) = self.evals {
            len_evals = evals.len();
        }
        len_coeffs.max(len_evals)
        // match self.coeffs {
        //     Some(ref coeffs) => coeffs.len(),
        //     None => match self.evals {
        //         Some(ref evals) => evals.len(),
        //         None => 0,
        //     },
        // }
    }

    fn __add__(&mut self, other: PolynomialOrGFElement) -> PyResult<Polynomial> {
        match other {
            PolynomialOrGFElement::Polynomial(mut b) => {
                // let mut a = self.clone();
                self.prepare_operation(&mut b)?;
                let evals: Vec<GFElement> = self
                    .evals
                    .as_ref()
                    .unwrap()
                    .par_iter()
                    .zip(b.evals.as_ref().unwrap())
                    .map(|(e1, e2)| e1.__add__(e2).unwrap())
                    .collect();
                Ok(Polynomial {
                    coeffs: None,
                    evals: Some(evals),
                    gf: self.gf.clone(),
                })
            }
            PolynomialOrGFElement::GFElement(b) => {
                // let mut a = self.clone();
                // self.calc_coeffs_if_necessary()?;
                let mut a = self.clone();
                // a.coeffs.as_mut().unwrap()[0] = a.coeffs.as_mut().unwrap()[0].__add__(&b)?;
                a.coeffs = match a.coeffs {
                    Some(ref mut coeffs) => {
                        coeffs[0] = coeffs[0].__add__(&b).unwrap();
                        Some(coeffs.to_vec())
                    }
                    None => None,
                };
                a.evals = match a.evals {
                    Some(ref evals) => {
                        Some(evals.par_iter().map(|e| e.__add__(&b).unwrap()).collect())
                    }
                    None => None,
                };
                Ok(a)
            }
        }
        // if isinstance(other, Polynomial):
        //     a, b = self._prepare_operation(self, other)
        //     assert a._evals is not None
        //     assert b._evals is not None
        //     return Polynomial(
        //         self._p,
        //         evals=[e1 + e2 for e1, e2 in zip(a._evals, b._evals)],
        //     )
        // elif isinstance(other, GFElement):
        //     assert self._p == other.gf.p
        //     if self._coeffs is not None:
        //         return Polynomial(
        //             self._p, coeffs=[self._coeffs[0] + other] + list(self._coeffs[1:])
        //         )
        //     elif self._evals is not None:
        //         return self + Polynomial(self._p, coeffs=[other])
        //     else:
        //         raise RuntimeError
        // else:
        //     raise RuntimeError
    }

    fn __sub__(&self, other: PolynomialOrGFElement) -> PyResult<Polynomial> {
        match other {
            PolynomialOrGFElement::Polynomial(mut b) => {
                let mut a = self.clone();
                a.prepare_operation(&mut b)?;
                let evals: Vec<GFElement> = a
                    .evals
                    .as_ref()
                    .unwrap()
                    .par_iter()
                    .zip(b.evals.unwrap())
                    .map(|(e1, e2)| e1.__sub__(&e2).unwrap())
                    .collect();
                Ok(Polynomial {
                    coeffs: None,
                    evals: Some(evals),
                    gf: a.gf,
                })
            }
            PolynomialOrGFElement::GFElement(b) => {
                let mut a = self.clone();
                // a.calc_coeffs_if_necessary()?;
                // a.coeffs.as_mut().unwrap()[0] = a.coeffs.as_mut().unwrap()[0].__sub__(&b)?;
                a.coeffs = match a.coeffs {
                    Some(ref mut coeffs) => {
                        coeffs[0] = coeffs[0].__sub__(&b).unwrap();
                        Some(coeffs.to_vec())
                    }
                    None => None,
                };
                a.evals = match a.evals {
                    Some(ref evals) => {
                        Some(evals.par_iter().map(|e| e.__sub__(&b).unwrap()).collect())
                    }
                    None => None,
                };
                Ok(a)
            }
        }
        // if isinstance(other, Polynomial):
        //     a, b = self._prepare_operation(self, other)
        //     assert a._evals is not None
        //     assert b._evals is not None
        //     return Polynomial(
        //         self._p,
        //         evals=[e1 - e2 for e1, e2 in zip(a._evals, b._evals)],
        //     )
        // elif isinstance(other, GFElement):
        //     assert self._p == other.gf.p
        //     if self._coeffs is not None:
        //         return Polynomial(
        //             self._p, coeffs=[self._coeffs[0] - other] + list(self._coeffs[1:])
        //         )
        //     elif self._evals is not None:
        //         return self - Polynomial(self._p, coeffs=[other])
        //     else:
        //         raise RuntimeError
        // else:
        //     raise RuntimeError
    }

    fn __mul__(&mut self, other: PolynomialOrGFElement) -> PyResult<Polynomial> {
        match other {
            PolynomialOrGFElement::Polynomial(mut b) => {
                // let mut a = self.clone();
                self.prepare_operation(&mut b)?;
                let evals: Vec<GFElement> = self
                    .evals
                    .as_ref()
                    .unwrap()
                    .par_iter()
                    .zip(b.evals.unwrap())
                    .map(|(e1, e2)| e1.__mul__(&e2).unwrap())
                    .collect();
                Ok(Polynomial {
                    coeffs: None,
                    evals: Some(evals),
                    gf: self.gf.clone(),
                })
            }
            PolynomialOrGFElement::GFElement(b) => {
                // self.calc_coeffs_if_necessary(None)?;
                // let mut a = self.clone();
                // a.coeffs = Some(
                //     a.coeffs
                //         .unwrap()
                //         .iter()
                //         .map(|e| e.__mul__(&b).unwrap())
                //         .collect(),
                // );
                // Ok(a)
                // self.calc_evals_if_necessary()?;
                let mut a = self.clone();
                // a.evals = Some(
                //     a.evals
                //         .unwrap()
                //         .iter()
                //         .map(|e| e.__mul__(&b).unwrap())
                //         .collect(),
                // );
                a.evals = match a.evals {
                    Some(ref evals) => {
                        Some(evals.par_iter().map(|e| e.__mul__(&b).unwrap()).collect())
                    }
                    None => None,
                };
                a.coeffs = match a.coeffs {
                    Some(ref coeffs) => {
                        Some(coeffs.par_iter().map(|e| e.__mul__(&b).unwrap()).collect())
                    }
                    None => None,
                };
                Ok(a)
            }
        }
        // if isinstance(other, Polynomial):
        //     a, b = self._prepare_operation(self, other)
        //     assert a._evals is not None
        //     assert b._evals is not None
        //     return Polynomial(
        //         self._p,
        //         evals=[e1 * e2 for e1, e2 in zip(a._evals, b._evals)],
        //     )
        // elif isinstance(other, GFElement):
        //     assert self._p == other.gf.p
        //     if self._coeffs is not None:
        //         return Polynomial(self._p, coeffs=[c * other for c in self._coeffs])
        //     elif self._evals is not None:
        //         return Polynomial(self._p, evals=[e * other for e in self._evals])
        //     else:
        //         raise RuntimeError
        // else:
        //     raise RuntimeError
    }

    fn __getitem__(&mut self, idx: usize) -> PyResult<GFElement> {
        // self._calc_coeffs_if_necessary()
        // assert self._coeffs is not None
        // return self._coeffs[idx]
        // let mut a = self.clone();
        self.calc_coeffs_if_necessary()?;
        if idx < self.__len__() {
            Ok(self.coeffs.as_ref().unwrap()[idx].clone())
        } else {
            Ok(self.gf.zero())
        }
    }

    fn degree(&mut self) -> PyResult<usize> {
        self.calc_coeffs_if_necessary()?;
        for (i, c) in self.coeffs.as_ref().unwrap().iter().enumerate().rev() {
            if !c.__eq__(self.gf.zero()) {
                return Ok(i);
            }
        }
        Ok(0)
        // self._calc_coeffs_if_necessary()
        // assert self._coeffs is not None
        // i = len(self._coeffs) - 1
        // while True:
        //     if self._coeffs[i] != self._Fp(0):
        //         return i
        //     i -= 1
        //     if i < 0:
        //         return 0
    }

    fn extend(&self, n: usize) -> PyResult<Polynomial> {
        let mut a = self.clone();
        a.calc_coeffs_if_necessary()?;
        if a.coeffs.as_ref().unwrap().len() >= n {
            return Ok(a);
        }
        for _ in 0..(n - a.coeffs.as_ref().unwrap().len()) {
            a.coeffs.as_mut().unwrap().push(self.gf.zero());
        }
        a.calc_evals_if_necessary()?;
        Ok(a)
        // assert n & (n - 1) == 0, "n must be power of 2"
        // self._calc_coeffs_if_necessary()
        // assert self._coeffs is not None
        // assert n > len(self._coeffs)
        // coeffs = list(self._coeffs) + [self._Fp(0)] * (n - len(self._coeffs))
        // self = Polynomial(self._p, coeffs=coeffs)
    }

    fn __divmod__(&mut self, other: &SparsePolynomial) -> PyResult<(Polynomial, Polynomial)> {
        let mut q = self.coeffs()?;
        let mut r = self.coeffs()?;
        let n = self.degree()?;
        let (max_idx, max_coeff) = &other.coeffs[other.coeffs.len() - 1];
        for idx in (*max_idx..n + 1).rev() {
            let k = r[idx].__truediv__(&max_coeff)?;
            q[idx - max_idx] = k.clone();
            for (j, c) in other.coeffs[0..other.coeffs.len() - 1].iter().rev() {
                r[idx - max_idx + j] = r[idx - max_idx + j].__sub__(&k.__mul__(&c)?)?;
            }
        }
        let q = q[0..n - max_idx + 1]
            .into_iter()
            .map(|x| BigIntOrGFElement::GFElement(x.clone()))
            .collect();
        let r = r[0..*max_idx]
            .into_iter()
            .map(|x| BigIntOrGFElement::GFElement(x.clone()))
            .collect();
        Ok((
            Polynomial::new(self.gf.clone(), Some(q), None)?,
            Polynomial::new(self.gf.clone(), Some(r), None)?,
        ))
    }

    #[getter]
    fn coeffs(&mut self) -> PyResult<Vec<GFElement>> {
        self.calc_coeffs_if_necessary()?;
        Ok(self.coeffs.as_ref().unwrap().clone())
    }

    #[getter]
    fn evals(&mut self) -> PyResult<Vec<GFElement>> {
        self.calc_evals_if_necessary()?;
        Ok(self.evals.as_ref().unwrap().clone())
    }
}

impl Polynomial {
    // fn fft(&self, coeffs: Vec<GFElement>, omegas: &Vec<GFElement>) -> PyResult<Vec<GFElement>> {
    //     let n = coeffs.len();
    //     if n & (n - 1) != 0 {
    //         return Err(PyValueError::new_err("n must be a power of 2"));
    //     }
    //     if n == 1 {
    //         return Ok(coeffs);
    //     }
    //     let ye = self.fft(coeffs.iter().step_by(2).cloned().collect(), omegas)?;
    //     let yo = self.fft(coeffs.iter().skip(1).step_by(2).cloned().collect(), omegas)?;
    //     let omega = &omegas[n.ilog2() as usize];
    //     let mut omega_i = self.gf.__call__(Value::GFElement(self.gf.one()));
    //     let mut evals = vec![self.gf.zero(); n];
    //     for i in 0..(n / 2) {
    //         let tmp = omega_i.__mul__(yo[i].clone())?;
    //         evals[i] = ye[i].__add__(tmp.clone())?;
    //         evals[i + n / 2] = ye[i].__sub__(tmp)?;
    //         omega_i = omega_i.__mul__(omega.clone())?;
    //     }
    //     Ok(evals)
    // }
    fn fft(&self, coeffs: &Vec<GFElement>, omegas: &[GFElement]) -> PyResult<Vec<GFElement>> {
        let n = coeffs.len();
        if n & (n - 1) != 0 {
            return Err(PyValueError::new_err("n must be a power of 2"));
        }
        // if n == 1 {
        //     return Ok(coeffs.clone());
        // }
        let k = n.ilog2();
        let step = MAX_L >> k;
        let mut evals = coeffs.clone();
        for (i, j) in REVBITS.iter().step_by(step as usize).enumerate() {
            if &(i as u32) < j {
                evals.swap(i, *j as usize);
            }
        }
        let mut r = 1;
        let one = self.gf.one();
        for omega in &omegas[1..(k + 1) as usize] {
            for l in (0..n).step_by(2 * r) {
                let mut omega_i = one.clone();
                for i in 0..r {
                    let tmp = omega_i.__mul__(&evals[l + i + r])?;
                    (evals[l + i], evals[l + i + r]) =
                        (evals[l + i].__add__(&tmp)?, evals[l + i].__sub__(&tmp)?);
                    omega_i = omega_i.__mul__(omega)?;
                }
            }
            r <<= 1;
        }
        Ok(evals)
    }

    fn coeffs_to_evals(&self, coeffs: &Vec<GFElement>) -> PyResult<Vec<GFElement>> {
        let mut omegas = Vec::new();
        let n = coeffs.len();
        let bit_length = n.ilog2() + 1;
        for i in 0..bit_length as usize {
            omegas.push(
                self.gf
                    .nth_root_of_unity(BigInt::from_str("2").unwrap().pow(i as u32))?,
            );
        }
        self.fft(coeffs, &omegas)
    }

    fn evals_to_coeffs(&self, evals: &Vec<GFElement>) -> PyResult<Vec<GFElement>> {
        let n = evals.len();
        if n & (n - 1) != 0 {
            return Err(PyValueError::new_err("n must be a power of 2"));
        }
        let ninv = self
            .gf
            .__call__(BigIntOrGFElement::BigInt(
                BigInt::from_str(&n.to_string()).unwrap(),
            ))
            .__pow__(BigInt::from_str("-1").unwrap(), None)?;
        let bit_length = n.ilog2() + 1;
        let mut omega_invs = Vec::new();
        for i in 0..bit_length as usize {
            omega_invs.push(
                (self
                    .gf
                    .nth_root_of_unity(BigInt::from_str("2").unwrap().pow(i as u32))?)
                .__pow__(BigInt::from_str("-1").unwrap(), None)?,
            );
        }
        let result = self.fft(evals, &omega_invs)?;
        Ok(result
            .par_iter()
            .map(|x| x.__mul__(&ninv).unwrap())
            .collect())
    }

    fn prepare_operation(&mut self, other: &mut Self) -> PyResult<()> {
        // TODO: __len__() がやばい。 evals と coeffs が違う値を持ってしまうとき (そもそもこんなことおきるわけなくない？誰が悪い？) に異常が起きる

        // println!("before");
        // if other.evals.is_none() {
        //     println!("{:?}", other.evals);
        // } else {
        //     println!("{:?}", other.evals.as_ref().unwrap().len());
        // }
        // self.calc_evals_if_necessary(Some(true))?;
        // other.calc_evals_if_necessary(Some(true))?;
        // self.calc_coeffs_if_necessary(Some(true))?;
        // other.calc_coeffs_if_necessary(Some(true))?;
        // println!("after");
        // if self.evals.is_none() {
        //     println!("{:?}", self.evals);
        // } else {
        //     println!("{:?}", self.evals.as_ref().unwrap().len());
        // }
        // if other.evals.is_none() {
        //     println!("{:?}", other.evals);
        // } else {
        //     println!("{:?}", other.evals.as_ref().unwrap().len());
        // }
        let len_self = self.__len__();
        let len_other = other.__len__();
        // println!(
        //     "before: self.len: {}, other.len: {}",
        //     self.__len__(),
        //     other.__len__()
        // );
        if len_self > len_other {
            other.extend_internal(len_self)?;
            // other.evals = Some(other.coeffs_to_evals(&other.coeffs.as_ref().unwrap())?);
            // other.coeffs = None;
        } else if len_self < len_other {
            self.extend_internal(len_other)?;
            // self.evals = Some(self.coeffs_to_evals(&self.coeffs.as_ref().unwrap())?);
            // self.coeffs = None;
        }
        self.calc_evals_if_necessary()?;
        other.calc_evals_if_necessary()?;
        if self.evals.as_ref().unwrap().len() != other.evals.as_ref().unwrap().len() {
            return Err(PyRuntimeError::new_err("Lengths of evals are different"));
        }
        // println!(
        //     "after: self.len: {}, other.len: {}",
        //     self.__len__(),
        //     other.__len__()
        // );
        Ok(())
    }

    fn extend_internal(&mut self, n: usize) -> PyResult<()> {
        self.calc_coeffs_if_necessary()?;
        if self.coeffs.as_ref().unwrap().len() >= n {
            return Ok(());
        }
        for _ in 0..(n - self.coeffs.as_ref().unwrap().len()) {
            self.coeffs.as_mut().unwrap().push(self.gf.zero());
        }
        if self.coeffs.as_ref().unwrap().len() != n {
            return Err(PyRuntimeError::new_err(
                "Lengths of coeffs are different in extend",
            ));
        }
        self.calc_evals_if_necessary()?;
        Ok(())
        // assert n & (n - 1) == 0, "n must be power of 2"
        // self._calc_coeffs_if_necessary()
        // assert self._coeffs is not None
        // assert n > len(self._coeffs)
        // coeffs = list(self._coeffs) + [self._Fp(0)] * (n - len(self._coeffs))
        // self = Polynomial(self._p, coeffs=coeffs)
    }

    fn calc_evals_if_necessary(&mut self) -> PyResult<()> {
        if self.evals.is_none()
            || (!self.coeffs.is_none()
                && self.evals.as_ref().unwrap().len() < self.coeffs.as_ref().unwrap().len())
        {
            match self.coeffs {
                Some(ref coeffs) => {
                    self.evals = Some(self.coeffs_to_evals(&coeffs)?);
                    // self.coeffs = None;
                }
                None => {
                    return Err(PyValueError::new_err("coeffs is None"));
                }
            }
        }
        Ok(())
    }

    fn calc_coeffs_if_necessary(&mut self) -> PyResult<()> {
        if self.coeffs.is_none()
            || (!self.evals.is_none()
                && self.coeffs.as_ref().unwrap().len() < self.evals.as_ref().unwrap().len())
        {
            match &self.evals {
                Some(evals) => {
                    self.coeffs = Some(self.evals_to_coeffs(&evals)?);
                    // self.evals = None;
                }
                None => {
                    return Err(PyValueError::new_err("coeffs is None"));
                }
            }
        }
        Ok(())
    }
}

impl<'a> FromPyObject<'a> for PolynomialOrGFElement {
    fn extract(ob: &'a PyAny) -> PyResult<Self> {
        if let Ok(s) = Polynomial::extract(ob) {
            Ok(PolynomialOrGFElement::Polynomial(s))
        } else if let Ok(s) = GFElement::extract(ob) {
            Ok(PolynomialOrGFElement::GFElement(s))
        } else {
            Err(PyValueError::new_err("Invalid type"))
        }
    }
}

#[pymethods]
impl SparsePolynomial {
    #[new]
    fn new(
        gf: GF,
        coeffs: Vec<(usize, BigIntOrGFElement)>,
        // evals: Option<Vec<BigIntOrGFElement>>,
    ) -> Self {
        let mut coeffs: Vec<(usize, GFElement)> = coeffs
            .iter()
            .map(|(i, c)| match c {
                BigIntOrGFElement::BigInt(c) => (*i, GFElement::new(c.clone(), gf.clone(), false)),
                BigIntOrGFElement::GFElement(c) => (*i, c.clone()),
            })
            .collect();
        coeffs.sort_by(|a, b| a.0.cmp(&b.0));
        Self { gf, coeffs }
    }

    fn __repr__(&self) -> PyResult<String> {
        let mut res = Vec::new();
        for (i, c) in &self.coeffs {
            if c.__eq__(self.gf.zero()) {
                continue;
            }
            if *i == 0 {
                res.push(c.__str__());
            } else if c.__eq__(self.gf.one()) {
                res.push(format!("x^{}", i));
            } else {
                res.push(format!("{} * x^{}", c.__str__(), i));
            }
        }
        let mut ret = res.join(" + ");
        if ret.len() > 1024 {
            ret = ret[..1024].to_string() + " ...";
        }
        Ok(format!("{} in {}", ret, self.gf.__str__()))
    }
}
