use std::{mem, ops::SubAssign};

use num_bigint::{BigInt, RandBigInt};
use num_traits::{One, Zero};
use rand;

#[derive(Clone, Debug)]
pub struct SS {
    /// threshold
    pub t: usize,
    /// total number of shares
    pub n: usize,
    /// prime in ff
    pub p: BigInt,
}

impl SS {
    pub fn split(&self, secret: BigInt) -> Vec<(usize, BigInt)> {
        assert!(self.t < self.n);
        let polynomial = self.sample_polynomial(secret);
        self.evaluate_polynomial(polynomial)
    }

    fn sample_polynomial(&self, secret: BigInt) -> Vec<BigInt> {
        let mut coeff: Vec<BigInt> = vec![secret];
        let mut rng = rand::thread_rng();
        let low = BigInt::zero();
        let high = &self.p - BigInt::one();
        let random_coeffs: Vec<BigInt> = (0..(self.t - 1))
            .map(|_| rng.gen_bigint_range(&low, &high))
            .collect();
        coeff.extend(random_coeffs);
        coeff
    }

    fn evaluate_polynomial(&self, polynomial: Vec<BigInt>) -> Vec<(usize, BigInt)> {
        (1..=self.n)
            .map(|x| (x, self.mod_evaluate_at(&polynomial, x)))
            .collect()
    }

    fn mod_evaluate_at(&self, ploynomial: &[BigInt], x: usize) -> BigInt {
        let x_bigint = BigInt::from(x);
        ploynomial
            .iter()
            .rev()
            .fold(Zero::zero(), |sum, item| (&x_bigint * sum + item) % &self.p)
    }

    pub fn recover(&self, shares: &[(usize, BigInt)]) -> BigInt {
        assert!(shares.len() == self.t, "wrong shares");
        let (xs, ys): (Vec<usize>, Vec<BigInt>) = shares.iter().cloned().unzip();
        let result = self.lagrange_interpolation(Zero::zero(), xs, ys);
        if result < Zero::zero() {
            result + &self.p
        } else {
            result
        }
    }

    fn lagrange_interpolation(&self, x: BigInt, xs: Vec<usize>, ys: Vec<BigInt>) -> BigInt {
        let len = xs.len();
        let xs_bigint: Vec<BigInt> = xs.iter().map(|x| BigInt::from(*x as i64)).collect();
        (0..len).fold(Zero::zero(), |sum, item| {
            let numerator = (0..len).fold(One::one(), |product: BigInt, i| {
                if i == item {
                    product
                } else {
                    product * (&x - &xs_bigint[i]) % &self.p
                }
            });
            let denominator = (0..len).fold(One::one(), |product: BigInt, i| {
                if i == item {
                    product
                } else {
                    product * (&xs_bigint[item] - &xs_bigint[i]) % &self.p
                }
            });
            (sum + numerator * self.mod_inv(denominator) * &ys[item]) % &self.p
        })
    }

    fn mod_inv(&self, a: BigInt) -> BigInt {
        let m = self.p.clone();
        let num = if a < Zero::zero() { a + &self.p } else { a };
        let (g, x, _) = SS::xgcd(num, m);
        assert!(g.is_one());
        (x + &self.p) % &self.p
    }

    fn xgcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
        let (mut sj, mut sj_last) = (BigInt::zero(), BigInt::one());
        let (mut tj, mut tj_last) = (BigInt::one(), BigInt::zero());
        let (mut rj, mut rj_last) = (b, a);

        while !rj.is_zero() {
            let quotient = &rj_last / &rj;
            rj_last.sub_assign(&quotient * &rj);
            sj_last.sub_assign(&quotient * &sj);
            tj_last.sub_assign(&quotient * &tj);
            mem::swap(&mut rj, &mut rj_last);
            mem::swap(&mut sj, &mut sj_last);
            mem::swap(&mut tj, &mut tj_last);
        }
        (rj_last, sj_last, tj_last)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recover_test() {
        let ss = SS {
            t: 3,
            n: 6,
            p: BigInt::from(1613),
        };

        let shares = ss.evaluate_polynomial(vec![
            BigInt::from(1234),
            BigInt::from(166),
            BigInt::from(94),
        ]);

        assert_eq!(
            shares,
            [
                (1, BigInt::from(1494)),
                (2, BigInt::from(329)),
                (3, BigInt::from(965)),
                (4, BigInt::from(176)),
                (5, BigInt::from(1188)),
                (6, BigInt::from(775)),
            ]
        );
        let r = ss.recover(&[
            (1, BigInt::from(1494)),
            (2, BigInt::from(329)),
            (3, BigInt::from(965)),
        ]);
        assert_eq!(r, BigInt::from(1234))
    }

    #[test]
    fn large_parime() {
        let ss = SS {
            t: 3,
            n: 6,
            p: BigInt::parse_bytes(
                b"fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f",
                16,
            )
            .unwrap(),
        };
        let secret = BigInt::parse_bytes(b"ffffffffffffffffffffffffffffffffffffff", 16).unwrap();
        let shares = ss.split(secret.clone());
        assert_eq!(secret, ss.recover(&shares[0..ss.t as usize]));
    }
}
