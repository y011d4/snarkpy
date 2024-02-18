use num_bigint::BigInt;
use once_cell::sync::Lazy;
use pyo3::{exceptions::PyValueError, prelude::*};
use std::str::FromStr;

use crate::field::{GFElement, Value, GF};

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
    #[pyo3(get)]
    coeffs: Option<Vec<GFElement>>,
    #[pyo3(get)]
    evals: Option<Vec<GFElement>>,
    gf: GF,
}

#[pymethods]
impl Polynomial {
    #[new]
    fn new(gf: GF, coeffs: Option<Vec<Value>>, evals: Option<Vec<Value>>) -> PyResult<Self> {
        if coeffs.is_none() && evals.is_none() {
            return Err(PyValueError::new_err(
                "coeffs and evals cannot both be None",
            ));
        }
        let coeffs = match coeffs {
            Some(coeffs) => {
                let mut result = Vec::new();
                for coeff in coeffs {
                    match coeff {
                        Value::BigInt(value) => {
                            result.push(GFElement::new(value, gf.clone(), false));
                        }
                        Value::GFElement(value) => {
                            result.push(value);
                        }
                    }
                }
                Some(result)
            }
            None => None,
        };
        let evals = match evals {
            Some(evals) => {
                let mut result = Vec::new();
                for eval in evals {
                    match eval {
                        Value::BigInt(value) => {
                            result.push(GFElement::new(value, gf.clone(), false));
                        }
                        Value::GFElement(value) => {
                            result.push(value);
                        }
                    }
                }
                Some(result)
            }
            None => None,
        };
        Ok(Polynomial { coeffs, evals, gf })
    }

    fn __call__(&mut self, x: Value) -> PyResult<GFElement> {
        self.calc_coeffs_if_necessary()?;
        let x = match x {
            Value::BigInt(value) => GFElement::new(value, self.gf.clone(), false),
            Value::GFElement(value) => value,
        };
        let mut result = self.gf.zero();
        for c in self.coeffs.as_ref().unwrap().iter().rev() {
            result = result
                .__mul__(x.clone())
                .unwrap()
                .__add__(c.clone())
                .unwrap();
        }
        Ok(result)
    }

    fn coeffs_to_evals(&mut self, coeffs: Vec<GFElement>) -> PyResult<Vec<GFElement>> {
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

    fn evals_to_coeffs(&mut self, evals: Vec<GFElement>) -> PyResult<Vec<GFElement>> {
        let n = evals.len();
        if n & (n - 1) != 0 {
            return Err(PyValueError::new_err("n must be a power of 2"));
        }
        let ninv = self
            .gf
            .__call__(Value::BigInt(BigInt::from_str(&n.to_string()).unwrap()))
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
            .iter()
            .map(|x| x.__mul__(ninv.clone()).unwrap())
            .collect())
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

    fn calc_evals_if_necessary(&mut self) -> PyResult<()> {
        if self.evals.is_none() {
            match self.coeffs {
                Some(ref coeffs) => {
                    self.evals = Some(self.coeffs_to_evals(coeffs.clone())?);
                }
                None => {
                    return Err(PyValueError::new_err("coeffs is None"));
                }
            }
        }
        Ok(())
    }

    fn calc_coeffs_if_necessary(&mut self) -> PyResult<()> {
        if self.coeffs.is_none() {
            match self.evals {
                Some(ref evals) => {
                    self.coeffs = Some(self.evals_to_coeffs(evals.clone())?);
                }
                None => {
                    return Err(PyValueError::new_err("coeffs is None"));
                }
            }
        }
        Ok(())
    }

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
        match self.coeffs {
            Some(ref coeffs) => coeffs.len(),
            None => match self.evals {
                Some(ref evals) => evals.len(),
                None => 0,
            },
        }
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
    fn fft(&self, coeffs: Vec<GFElement>, omegas: &Vec<GFElement>) -> PyResult<Vec<GFElement>> {
        let n = coeffs.len();
        if n & (n - 1) != 0 {
            return Err(PyValueError::new_err("n must be a power of 2"));
        }
        if n == 1 {
            return Ok(coeffs);
        }
        let k = n.ilog2();
        let step = MAX_L >> k;
        let mut evals = coeffs.clone();
        for (i, j) in REVBITS.iter().step_by(step as usize).enumerate() {
            evals.swap(i, *j as usize);
        }
        let mut r = 1;
        for omega in &omegas[1..(k + 1) as usize] {
            for l in (0..n).step_by(2 * r) {
                let mut omega_i = self.gf.__call__(Value::GFElement(self.gf.one()));
                for i in 0..r {
                    let tmp = omega_i.__mul__(coeffs[l + i + r].clone())?;
                    (evals[l + i], evals[l + i + r]) = (
                        coeffs[l + i].__add__(tmp.clone())?,
                        coeffs[l + i].__sub__(tmp)?,
                    );
                    omega_i = omega_i.__mul__(omega.clone())?;
                }
            }
            r <<= 1;
        }
        Ok(evals)
    }
}
