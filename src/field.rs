use num_bigint::BigUint;
use num_traits::{One, Zero};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldElement {
    pub num: BigUint,
    pub prime: BigUint,
}

impl FieldElement {
    pub fn new(num: BigUint, prime: BigUint) -> Result<Self, &'static str> {
        if prime.is_zero() {
            return Err("Prime cannot be zero");
        }
        let num_mod_prime = num % &prime;
        Ok(Self {
            num: num_mod_prime,
            prime,
        })
    }

    // Check if they are on the same prime field.
    fn check_same_prime(&self, other: &Self) -> Result<(), &'static str> {
        if self.prime != other.prime {
            Err("Different prime numbers")
        } else {
            Ok(())
        }
    }

    // Add two field elements.
    pub fn add(&self, other: &Self) -> Result<Self, &'static str> {
        self.check_same_prime(other)?;
        let sum = (&self.num + &other.num) % &self.prime;
        Self::new(sum, self.prime.clone())
    }

    /// Subtract two field elements.
    pub fn sub(&self, other: &Self) -> Result<Self, &'static str> {
        self.check_same_prime(other)?;
        let diff = (&self.num + &self.prime - &other.num) % &self.prime;
        Self::new(diff, self.prime.clone())
    }

    /// Multiply two field elements.
    pub fn mul(&self, other: &Self) -> Result<Self, &'static str> {
        self.check_same_prime(other)?;
        let product = (&self.num * &other.num) % &self.prime;
        Self::new(product, self.prime.clone())
    }

    /// Divide two field elements.
    pub fn div(&self, other: &Self) -> Result<Self, &'static str> {
        let inv_other = other.inv()?;
        self.mul(&inv_other)
    }

    /// Calculate the power of a field element.
    pub fn pow(&self, exp: u32) -> Result<Self, &'static str> {
        let result = self.num.modpow(&BigUint::from(exp), &self.prime);
        Self::new(result, self.prime.clone())
    }

    /// Calculate the inverse of a field element.
    pub fn inv(&self) -> Result<Self, &'static str> {
        if self.is_zero() {
            return Err("Cannot compute inverse of zero");
        }
        // Use Fermat's little theorem
        let inv = self.num.modpow(&(self.prime.clone() - BigUint::from(2u32)), &self.prime);
        Self::new(inv, self.prime.clone())
    }

    /// Calculate the additive inverse of a field element.
    pub fn negate(&self) -> Result<Self, &'static str> {
        Self::zero(&self.prime).sub(self)
    }

    /// Returns a zero element.
    pub fn zero(prime: &BigUint) -> Self {
        Self {
            num: BigUint::zero(),
            prime: prime.clone(),
        }
    }

    /// Returns a one element.
    pub fn one(prime: &BigUint) -> Self {
        Self {
            num: BigUint::one(),
            prime: prime.clone(),
        }
    }

    /// Check if the element is zero.
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    #[test]
    fn test_add() {
        let prime = BigUint::from(7u32);
        let a = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(5u32), prime.clone()).unwrap();
        let result = a.add(&b).unwrap();
        assert_eq!(result.num, BigUint::from(1u32));
    }

    #[test]
    fn test_sub() {
        let prime = BigUint::from(7u32);
        let a = FieldElement::new(BigUint::from(2u32), prime.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(5u32), prime.clone()).unwrap();
        let result = a.sub(&b).unwrap();
        assert_eq!(result.num, BigUint::from(4u32));
    }

    #[test]
    fn test_mul() {
        let prime = BigUint::from(7u32);
        let a = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(5u32), prime.clone()).unwrap();
        let result = a.mul(&b).unwrap();
        assert_eq!(result.num, BigUint::from(1u32));
    }
    
    #[test]
    fn test_inv() {
        let prime = BigUint::from(7u32);
        let a = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        let inv_a = a.inv().unwrap();
        assert_eq!(inv_a.num, BigUint::from(5u32));
    }

    #[test]
    fn test_div() {
        let prime = BigUint::from(7u32);
        let a = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(5u32), prime.clone()).unwrap();
        let result = a.div(&b).unwrap();
        assert_eq!(result.num, BigUint::from(2u32));
    }

    #[test]
    fn test_negate() {
        let prime = BigUint::from(7u32);
        let a = FieldElement::new(BigUint::from(2u32), prime.clone()).unwrap();
        let neg_a = a.negate().unwrap();
        assert_eq!(neg_a.num, BigUint::from(5u32));
    }
}
