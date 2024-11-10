use num_bigint::BigUint;
use num_traits::{One, Zero};
use num_traits::{FromPrimitive, ToPrimitive};

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

    pub fn find_primitive_root(order: usize, prime: &BigUint) -> Result<BigUint, &'static str> {
        use num_traits::One;
        let one = BigUint::one();
        let p_minus_one = prime - &one;
        let order_biguint = BigUint::from(order);
    
        if &p_minus_one % &order_biguint != BigUint::zero() {
            return Err("Order does not divide p - 1");
        }
    
        let cofactor = &p_minus_one / &order_biguint;
    
        let mut factors = Vec::new();
        let mut n = order;
        let mut i = 2;
        while i * i <= n {
            if n % i == 0 {
                factors.push(i);
                while n % i == 0 {
                    n /= i;
                }
            }
            i += 1;
        }
        if n > 1 {
            factors.push(n);
        }
    
        for candidate in 2..prime.to_usize().unwrap() {
            let candidate_biguint = BigUint::from(candidate);
    
            if candidate_biguint.modpow(&cofactor, prime) == one {
                continue;
            }
    
            let mut is_generator = true;
            for &factor in &factors {
                let exp = &p_minus_one / BigUint::from(factor);
                if candidate_biguint.modpow(&exp, prime) == one {
                    is_generator = false;
                    break;
                }
            }
    
            if is_generator {
                return Ok(candidate_biguint.modpow(&cofactor, prime));
            }
        }
    
        Err("No primitive root found")
    }

    pub fn from_u32(n: u32, prime: &BigUint) -> Result<Self, &'static str> {
        let num = BigUint::from_u32(n).ok_or("Invalid u32")?;
        FieldElement::new(num, prime.clone())
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

    #[test]
    fn test_find_primitive_root_and_subgroup() {
        let prime = BigUint::from(31u32);
        let order = 5;

        let primitive_root = FieldElement::find_primitive_root(order, &prime).unwrap();
        let mut group = Vec::new();

        for i in 0..order {
            let power = primitive_root.modpow(&BigUint::from(i), &prime);
            let element = FieldElement::new(power, prime.clone()).unwrap();
            group.push(element);
        }

        let expected = vec![
            FieldElement::new(BigUint::from(1u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(2u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(4u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(8u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(16u32), prime.clone()).unwrap(),
        ];

        assert_eq!(group, expected);
    }
}
