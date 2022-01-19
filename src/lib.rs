use num_bigint::{BigInt, RandBigInt};
use num_traits::{One, Zero};
use rand;

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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
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
        )
    }
}
