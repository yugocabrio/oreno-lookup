use std::collections::HashMap;
use crate::field::FieldElement;
use num_bigint::BigUint;

/// 多変数多項式を表現する構造体
#[derive(Clone, Debug, PartialEq)]
pub struct MultivariatePolynomial {
    /// 単項のマップ: (変数の指数ベクトル) => 係数
    pub terms: HashMap<Vec<u32>, FieldElement>,
}

impl MultivariatePolynomial {
    /// 新しい多変数多項式を作成
    pub fn new(terms: HashMap<Vec<u32>, FieldElement>) -> Self {
        Self { terms }
    }

    /// 多変数多項式の加算
    pub fn multi_poly_add(&self, other: &Self) -> Result<Self, &'static str> {
        let prime = self.terms.values().next().unwrap().prime.clone();
        let mut result_terms = self.terms.clone();

        for (exponent, coeff) in &other.terms {
            let entry = result_terms.entry(exponent.clone()).or_insert(FieldElement::zero(&prime));
            *entry = entry.add(coeff)?;
        }

        Ok(Self::new(result_terms))
    }

    /// 多変数多項式の減算
    pub fn multi_poly_sub(&self, other: &Self) -> Result<Self, &'static str> {
        let prime = self.terms.values().next().unwrap().prime.clone();
        let mut result_terms = self.terms.clone();

        for (exponent, coeff) in &other.terms {
            let entry = result_terms.entry(exponent.clone()).or_insert(FieldElement::zero(&prime));
            *entry = entry.sub(coeff)?;
        }

        Ok(Self::new(result_terms))
    }

    /// 多変数多項式の乗算
    pub fn multi_poly_mul(&self, other: &Self) -> Result<Self, &'static str> {
        let prime = self.terms.values().next().unwrap().prime.clone();
        let mut result_terms = HashMap::new();

        for (exp1, coeff1) in &self.terms {
            for (exp2, coeff2) in &other.terms {
                // 指数ベクトルの加算
                let mut new_exp = vec![];
                let len = usize::max(exp1.len(), exp2.len());
                for i in 0..len {
                    let e1 = if i < exp1.len() { exp1[i] } else { 0 };
                    let e2 = if i < exp2.len() { exp2[i] } else { 0 };
                    new_exp.push(e1 + e2);
                }
                let coeff_product = coeff1.mul(coeff2)?;
                let entry = result_terms.entry(new_exp).or_insert(FieldElement::zero(&prime));
                *entry = entry.add(&coeff_product)?;
            }
        }

        Ok(Self::new(result_terms))
    }

    /// 多変数多項式の評価
    pub fn multi_poly_evaluate(&self, values: &[FieldElement]) -> Result<FieldElement, &'static str> {
        let prime = self.terms.values().next().unwrap().prime.clone();
        let mut result = FieldElement::zero(&prime);

        for (exponents, coeff) in &self.terms {
            let mut term_value = coeff.clone();
            for (i, &exp) in exponents.iter().enumerate() {
                if i >= values.len() {
                    return Err("Not enough values provided for evaluation");
                }
                let value_pow_exp = values[i].pow(exp)?; // 修正ポイント
                term_value = term_value.mul(&value_pow_exp)?;
            }
            result = result.add(&term_value)?;
        }

        Ok(result)
    }

    /// 多変数多項式の変数に関する微分
    pub fn muti_poly_derivative(&self, var_index: usize) -> Result<Self, &'static str> {
        let prime = self.terms.values().next().unwrap().prime.clone();
        let mut result_terms = HashMap::new();

        for (exponents, coeff) in &self.terms {
            if var_index >= exponents.len() || exponents[var_index] == 0 {
                continue;
            }
            let mut new_exponents = exponents.clone();
            let exp_k = new_exponents[var_index];
            new_exponents[var_index] -= 1;
            let exp_k_fe = FieldElement::from_u32(exp_k, &prime)?; // 修正ポイント
            let new_coeff = coeff.mul(&exp_k_fe)?;
            let entry = result_terms.entry(new_exponents).or_insert(FieldElement::zero(&prime));
            *entry = entry.add(&new_coeff)?;
        }

        Ok(Self::new(result_terms))
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::FieldElement;
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;

    fn field_element(n: u32, prime: u32) -> FieldElement {
        FieldElement::from_u32(n, &BigUint::from_u32(prime).unwrap()).unwrap()
    }

    #[test]
    fn test_multivariate_polynomial_add() {
        let prime = 17;
        let coeff1 = field_element(3, prime);
        let coeff2 = field_element(5, prime);

        let mut terms1 = HashMap::new();
        terms1.insert(vec![1, 0], coeff1.clone());

        let mut terms2 = HashMap::new();
        terms2.insert(vec![1, 0], coeff2.clone());

        let poly1 = MultivariatePolynomial::new(terms1);
        println!("poly1: {:?}", poly1);
        let poly2 = MultivariatePolynomial::new(terms2);
        println!("poly2: {:?}", poly2);

        let result = poly1.multi_poly_add(&poly2).unwrap();
        println!("result: {:?}", result);

        let mut expected_terms = HashMap::new();
        expected_terms.insert(vec![1, 0], coeff1.add(&coeff2).unwrap());

        let expected_poly = MultivariatePolynomial::new(expected_terms);

        assert_eq!(result, expected_poly);
    }

    #[test]
    fn test_multivariate_polynomial_mul() {
        let prime = 17;
        let coeff1 = field_element(2, prime);
        let coeff2 = field_element(3, prime);

        let mut terms1 = HashMap::new();
        terms1.insert(vec![1, 0], coeff1.clone()); // 2 * x^1 * y^0

        let mut terms2 = HashMap::new();
        terms2.insert(vec![0, 1], coeff2.clone()); // 3 * x^0 * y^1

        let poly1 = MultivariatePolynomial::new(terms1);
        println!("poly1: {:?}", poly1);
        let poly2 = MultivariatePolynomial::new(terms2);
        println!("poly2: {:?}", poly2);

        let result = poly1.multi_poly_mul(&poly2).unwrap();
        println!("poly3: {:?}", result);

        let mut expected_terms = HashMap::new();
        expected_terms.insert(vec![1, 1], coeff1.mul(&coeff2).unwrap()); // 6 * x^1 * y^1

        let expected_poly = MultivariatePolynomial::new(expected_terms);

        assert_eq!(result, expected_poly);
    }

    #[test]
    fn test_multivariate_polynomial_evaluate() {
        let prime = 17;
        let coeff = field_element(4, prime);

        let mut terms = HashMap::new();
        terms.insert(vec![1, 2], coeff.clone()); // 4 * x^1 * y^2

        let poly = MultivariatePolynomial::new(terms);

        let x = field_element(2, prime);
        let y = field_element(3, prime);

        let result = poly.multi_poly_evaluate(&[x, y]).unwrap();

        // 計算: 4 * (2)^1 * (3)^2 = 4 * 2 * 9 = 72 mod 17 = 4
        let expected = field_element(4, prime);

        assert_eq!(result, expected);
    }

    #[test]
    fn test_multivariate_polynomial_derivative() {
        let prime = 17;
        let coeff = field_element(5, prime);

        let mut terms = HashMap::new();
        terms.insert(vec![2, 1], coeff.clone()); // 5 * x^2 * y^1

        let poly = MultivariatePolynomial::new(terms);

        // x に関する微分
        let derivative = poly.muti_poly_derivative(0).unwrap();

        let expected_coeff = coeff.mul(&field_element(2, prime)).unwrap(); // 5 * 2 = 10

        let mut expected_terms = HashMap::new();
        expected_terms.insert(vec![1, 1], expected_coeff); // 10 * x^1 * y^1

        let expected_poly = MultivariatePolynomial::new(expected_terms);

        assert_eq!(derivative, expected_poly);
    }
}
