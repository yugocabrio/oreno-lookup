use crate::field::FieldElement;

#[derive(Debug, Clone, PartialEq)]
pub struct UnivariatePolynomial {
    pub coefficients: Vec<FieldElement>,
}

impl From<Vec<FieldElement>> for UnivariatePolynomial {
    fn from(coefficients: Vec<FieldElement>) -> Self {
        Self::new(coefficients)
    }
}

impl UnivariatePolynomial {
    /// Creates a new univariate polynomial.
    pub fn new(coefficients: Vec<FieldElement>) -> Self {
        Self { coefficients }
    }

    /// Adds two univariate polynomials.
    pub fn uni_poly_add(&self, other: &Self) -> Result<Self, &'static str> {
        let len = usize::max(self.coefficients.len(), other.coefficients.len());
        let mut result_coeffs = Vec::with_capacity(len);
        let prime = self.coefficients[0].prime.clone();

        for i in 0..len {
            let a = self.coefficients.get(i).cloned().unwrap_or_else(|| FieldElement::zero(&prime));
            let b = other.coefficients.get(i).cloned().unwrap_or_else(|| FieldElement::zero(&prime));
            result_coeffs.push(a.add(&b)?);
        }

        Ok(Self::new(result_coeffs))
    }

    /// Subtracts another univariate polynomial from one univariate polynomial.
    pub fn uni_poly_sub(&self, other: &Self) -> Result<Self, &'static str> {
        let len = usize::max(self.coefficients.len(), other.coefficients.len());
        let mut result_coeffs = Vec::with_capacity(len);
        let prime = self.coefficients[0].prime.clone();

        for i in 0..len {
            let a = self.coefficients.get(i).cloned().unwrap_or_else(|| FieldElement::zero(&prime));
            let b = other.coefficients.get(i).cloned().unwrap_or_else(|| FieldElement::zero(&prime));
            result_coeffs.push(a.sub(&b)?);
        }

        Ok(Self::new(result_coeffs))
    }

    /// Multiplies two univariate polynomials.
    pub fn uni_poly_mul(&self, other: &Self) -> Result<Self, &'static str> {
        let len = self.coefficients.len() + other.coefficients.len() - 1;
        let mut result_coeffs = vec![FieldElement::zero(&self.coefficients[0].prime); len];

        for (i, a_coeff) in self.coefficients.iter().enumerate() {
            for (j, b_coeff) in other.coefficients.iter().enumerate() {
                let product = a_coeff.mul(b_coeff)?;
                result_coeffs[i + j] = result_coeffs[i + j].add(&product)?;
            }
        }

        Ok(Self::new(result_coeffs))
    }

    /// Scales the univariate polynomial by a scalar.
    pub fn uni_poly_scale(&self, scalar: &FieldElement) -> Result<Self, &'static str> {
        let scaled_coeffs = self
            .coefficients
            .iter()
            .map(|c| c.mul(scalar))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self::new(scaled_coeffs))
    }

    /// Evaluates the univariate polynomial at a given point.
    pub fn uni_poly_evaluate(&self, x: &FieldElement) -> Result<FieldElement, &'static str> {
        let mut result = FieldElement::zero(&x.prime);
        let mut x_pow = FieldElement::one(&x.prime);
        for coeff in &self.coefficients {
            result = result.add(&coeff.mul(&x_pow)?)?;
            // updte the power of x.
            x_pow = x_pow.mul(x)?;
        }
        Ok(result)
    }
    

    /// Composes one univariate polynomial with another univariate polynomial.
    pub fn uni_poly_compose(&self, other: &Self) -> Result<Self, &'static str> {
        let mut result = UnivariatePolynomial::new(vec![FieldElement::zero(&self.coefficients[0].prime)]);

        for coeff in self.coefficients.iter().rev() {
            result = result.uni_poly_mul(other)?;
            result = result.uni_poly_add(&UnivariatePolynomial::new(vec![coeff.clone()]))?;
        }
        result.trim();
        Ok(result)
    }

    /// Trims the univariate polynomial by removing leading zero coefficients.
    pub fn trim(&mut self) {
        while self.coefficients.last().map_or(false, |c| c.is_zero()) {
            self.coefficients.pop();
        }
    }

    /// Raises the univariate polynomial to a given power.
    pub fn poly_pow(&self, exponent: u32) -> Result<Self, &'static str> {
        let mut result = UnivariatePolynomial::new(vec![FieldElement::one(&self.coefficients[0].prime)]);
        let mut base = self.clone();
        let mut exp = exponent;

        while exp > 0 {
            if exp % 2 == 1 {
                result = result.uni_poly_mul(&base)?;
            }
            base = base.uni_poly_mul(&base)?;
            exp /= 2;
        }

        Ok(result)
    }

    /// Constructs an univariate polynomial using Lagrange interpolation.
    /// The polynomial is expressed as p(x) = Σ(y_i * L_i(x)) for i = 0 to n-1, where:
    /// L_i(x) = Π((x - x_j) / (x_i - x_j)) for k ≠ j, representing the Lagrange basis polynomials.
    pub fn uni_poly_lagrange_interpolation(
        // domain
        x_values: &[FieldElement],
        // evaluation
        y_values: &[FieldElement],
    ) -> Result<Self, &'static str> {
        if x_values.len() != y_values.len() {
            return Err("x_values and y_values must have the same length");
        }

        let prime = x_values[0].prime.clone();
        let mut result = UnivariatePolynomial::new(vec![FieldElement::zero(&prime)]);

        for (i, _x_i) in x_values.iter().enumerate() {
            let y_i = y_values[i].clone();
            // lagrange basis polynomia; L_i(x)
            let lagrange_basis = UnivariatePolynomial::uni_poly_lagrange_basis(x_values, i, &_x_i)?;
            // y_i * L_i(x)
            let term = lagrange_basis.uni_poly_scale(&y_i)?;
            result = result.uni_poly_add(&term)?;
        }

        Ok(result)
    }

    /// L_i(x) = Π((x - x_j) / (x_i - x_j)) for j ≠ i, representing the Lagrange basis polynomials.
    pub fn uni_poly_lagrange_basis(
        H: &[FieldElement],
        index: usize,
        _n_field: &FieldElement,
    ) -> Result<UnivariatePolynomial, &'static str> {
        let prime = H[0].prime.clone();
        let x_i = &H[index];

        // (x - x_j)
        let mut numerator = UnivariatePolynomial::new(vec![FieldElement::one(&prime)]);
        // (x_i - x_j)
        let mut denominator = FieldElement::one(&prime);

        for (j, x_j) in H.iter().enumerate() {
            if j != index {
                // (x - x_j)
                let term = UnivariatePolynomial::new(vec![x_j.negate()?, FieldElement::one(&prime)]);
                numerator = numerator.uni_poly_mul(&term)?;
                // (x_i - x_j)
                let denom_term = x_i.sub(x_j)?;
                denominator = denominator.mul(&denom_term)?;
            }
        }

        // (denominator)^{-1}
        let denominator_inv = denominator.inv()?;
        // L_i(x) = numerator / denominator
        numerator.uni_poly_scale(&denominator_inv)
    }


    /// Checks if the univariate polynomial is the zero polynomial.
    pub fn uni_poly_is_zero(&self) -> bool {
        self.coefficients.iter().all(|c| c.is_zero())
    }

    /// Returns the degree of the univariate polynomial.
    pub fn uni_poly_degree(&self) -> usize {
        let mut degree = self.coefficients.len();
        while degree > 0 && self.coefficients[degree - 1].is_zero() {
            degree -= 1;
        }
        degree.saturating_sub(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::FieldElement;
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use num_traits::{One, Zero};

    fn field_element(n: u32, prime: u32) -> FieldElement {
        FieldElement::new(
            BigUint::from_u32(n).unwrap(),
            BigUint::from_u32(prime).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn test_univariate_polynomial_add() {
        // 1 + 2x
        let a = UnivariatePolynomial::new(vec![
            field_element(1, 17),
            field_element(2, 17),
        ]);
        // 3 + 4x + 5x^2
        let b = UnivariatePolynomial::new(vec![
            field_element(3, 17),
            field_element(4, 17),
            field_element(5, 17),
        ]);
        // 4 + 6x + 5x^2
        let result = a.uni_poly_add(&b).unwrap();
        let expected = UnivariatePolynomial::new(vec![
            field_element(4, 17),
            field_element(6, 17),
            field_element(5, 17),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_univariate_polynomial_sub() {
        // 5 + 7x
        let a = UnivariatePolynomial::new(vec![
            field_element(5, 17),
            field_element(7, 17),
        ]);
        // 3 + 4x + x^2
        let b = UnivariatePolynomial::new(vec![
            field_element(3, 17),
            field_element(4, 17),
            field_element(1, 17),
        ]); 
        let result = a.uni_poly_sub(&b).unwrap();
        // 2 + 3x + 16x^2
        let expected = UnivariatePolynomial::new(vec![
            field_element(2, 17),
            field_element(3, 17),
            field_element(16, 17),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_univariate_polynomial_mul() {
        // 1 + 2x
        let a = UnivariatePolynomial::new(vec![
            field_element(1, 17),
            field_element(2, 17),
        ]);
        // 3 + 4x
        let b = UnivariatePolynomial::new(vec![
            field_element(3, 17),
            field_element(4, 17),
        ]);
        let result = a.uni_poly_mul(&b).unwrap();
        // 3 + 10x^2 + 8x^3
        let expected = UnivariatePolynomial::new(vec![
            field_element(3, 17),
            field_element(10, 17),
            field_element(8, 17),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_univariate_polynomial_scale() {
        // 2 + 3x
        let poly = UnivariatePolynomial::new(vec![
            field_element(2, 17),
            field_element(3, 17),
        ]);
        // 4
        let scalar = field_element(4, 17);
        let result = poly.uni_poly_scale(&scalar).unwrap();
        // 8 + 12x
        let expected = UnivariatePolynomial::new(vec![
            field_element(8, 17),
            field_element(12, 17),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_univariate_polynomial_evaluate() {
        // 2 + 3x + x^2
        let poly = UnivariatePolynomial::new(vec![
            field_element(2, 17),
            field_element(3, 17),
            field_element(1, 17),
        ]);

        // 4
        let x = field_element(4, 17);
        let result = poly.uni_poly_evaluate(&x).unwrap();

        // 30 mod 17 = 13
        let expected = field_element(13, 17);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_univariate_polynomial_trim() {
        // 2 + 3x + 0x^2 + 0x^3
        let mut poly = UnivariatePolynomial::new(vec![
            field_element(2, 17),
            field_element(3, 17),
            field_element(0, 17),
            field_element(0, 17),
        ]);
        poly.trim();

        let expected = UnivariatePolynomial::new(vec![
            field_element(2, 17),
            field_element(3, 17),
        ]);
        assert_eq!(poly, expected);
    }

    #[test]
    fn test_univariate_polynomial_compose() {
        let prime = BigUint::from_u32(17u32).unwrap();

        // f(x) = x + 1
        let one = FieldElement::new(BigUint::one(), prime.clone()).unwrap();
        let f = UnivariatePolynomial::new(vec![one.clone(), one.clone()]);

        // g(x) = x^2
        let zero = FieldElement::new(BigUint::zero(), prime.clone()).unwrap();
        let g = UnivariatePolynomial::new(vec![zero.clone(), zero.clone(), one.clone()]);

        // f(g(x))
        let composed = f.uni_poly_compose(&g).unwrap();

        // x^2 + 1
        let expected = UnivariatePolynomial::new(vec![one.clone(), zero.clone(), one.clone()]);
        assert_eq!(composed, expected);
    }

    #[test]
    fn test_univariate_polynomial_pow() {
        // 1 + x
        let poly = UnivariatePolynomial::new(vec![
            field_element(1, 17),
            field_element(1, 17),
        ]);
        // (1 + x)^3
        let powered = poly.poly_pow(3).unwrap();
        // 1 + 3x + 3x^2 + x^3
        let expected = UnivariatePolynomial::new(vec![
            field_element(1, 17),
            field_element(3, 17),
            field_element(3, 17),
            field_element(1, 17), 
        ]);
        assert_eq!(powered, expected);
    }

    #[test]
    // https://blog.lambdaclass.com/diving-deep-fri/
    fn test_univariate_polynomial_lagrange_interpolation() {
        let x_values = vec![
            field_element(1, 17),
            field_element(13, 17),
            field_element(16, 17),
            field_element(4, 17),
        ];
        let y_values = vec![
            field_element(3, 17),
            field_element(9, 17),
            field_element(13, 17),
            field_element(16, 17),
        ];

        let poly = UnivariatePolynomial::uni_poly_lagrange_interpolation(&x_values, &y_values).unwrap();

        let expected_poly = vec![
            field_element(6, 17),
            field_element(16, 17),
            field_element(2, 17),
            field_element(13, 17),
        ];

        assert_eq!(poly.coefficients, expected_poly);
    }

    #[test]
    fn test_univariate_polynomial_degree() {
        // 2 + 3x + 0x^2
        let poly = UnivariatePolynomial::new(vec![
            field_element(2, 17),
            field_element(3, 17),
            field_element(0, 17),
        ]);
        assert_eq!(poly.uni_poly_degree(), 1);
    }

    #[test]
    fn test_univariate_compute_lagrange_basis() {
        // p = 17
        let prime = BigUint::from_u32(17u32).unwrap();

        // generator g = 3
        let g = field_element(3, 17);
        // cyclic Subgroup H
        let mut H = Vec::new();
        // order 4
        let order = 4;

        for i in 1..=order {
            let g_i = g.pow(i as u32).unwrap();
            H.push(g_i);
        }

        println!("Cyclic Subgroup H: {:?}", H);

        // Lagrange basis polynomials
        for (i, x_i) in H.iter().enumerate() {
            let index = i;
            let L_i = UnivariatePolynomial::uni_poly_lagrange_basis(&H, index, &x_i);

            println!("L_{}: Lagrange basis polynomials: {:?}", i + 1, L_i);

            // evaluate L_i at one's in H
            for (j, x_j) in H.iter().enumerate() {
                let evaluated = L_i.as_ref().unwrap().uni_poly_evaluate(x_j).expect("Polynomial evaluation failed");
                if i == j {
                    // L_i(g^i) = 1
                    assert_eq!(evaluated, FieldElement::one(&prime));
                } else {
                    // L_i(g^j) = 0 (j ≠ i)
                    assert_eq!(evaluated, FieldElement::zero(&prime));
                }
            }
        }
    }
}
