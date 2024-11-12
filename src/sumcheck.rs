use crate::field::FieldElement;
use crate::multivariate_polynomial::MultivariatePolynomial;
use crate::univariate_polynomial::UnivariatePolynomial;
use num_bigint::BigUint;

pub struct Sumcheck<'a> {
    pub polynomial: &'a MultivariatePolynomial,
    pub num_variables: usize,
    pub verifier_challenges: Vec<FieldElement>,
}

impl<'a> Sumcheck<'a> {
    pub fn new(
        polynomial: &'a MultivariatePolynomial,
        num_variables: usize,
        verifier_challenges: Vec<FieldElement>,
    ) -> Self {
        Self {
            polynomial,
            num_variables,
            verifier_challenges,
        }
    }

    pub fn execute(&self) -> Result<(), &'static str> {
        let prime = &self.polynomial.terms.values().next().unwrap().prime;
        let mut claimed_sum = self.compute_initial_sum()?;
        println!("Initial total sum (H): {:?}", claimed_sum);
        let mut previous_r_i = Vec::new();

        for i in 0..self.num_variables {
            // Prover computes univariate polynomial
            let univariate_poly = self.compute_univariate_polynomial(i, &previous_r_i)?;
            println!(
                "Round {}: univariate polynomial g{}(x{}) = {:?}",
                i + 1,
                i + 1,
                i + 1,
                univariate_poly
            );

            // Evaluations at 0 and 1
            let x0 = FieldElement::zero(prime);
            let x1 = FieldElement::one(prime);
            let s_0 = univariate_poly.uni_poly_evaluate(&x0)?;
            let s_1 = univariate_poly.uni_poly_evaluate(&x1)?;
            println!(
                "Evaluation of univariate polynomial: S_0 = {:?}, S_1 = {:?}",
                s_0, s_1
            );

            // Verifier checks that s_0 + s_1 == claimed_sum
            let s_total = s_0.add(&s_1)?;
            if claimed_sum != s_total {
                return Err("Consistency check failed in Sumcheck protocol");
            }
            println!(
                "S_0 + S_1 = Claimed Sum = {:?}",
                claimed_sum
            );

            // Verifier provides challenge
            let r_i = self
                .verifier_challenges
                .get(i)
                .ok_or("Not enough challenges provided")?;
            println!("Random challenge r{} = {:?}", i + 1, r_i);

            // Prover and Verifier update claimed sum
            claimed_sum = univariate_poly.uni_poly_evaluate(r_i)?;
            println!("Updated Sum：{:?}", claimed_sum);

            previous_r_i.push(r_i.clone());
        }

        // Final check
        let final_value = self.polynomial.multi_poly_evaluate(&previous_r_i)?;
        println!("Final evaluation：{:?}", final_value);
        if claimed_sum != final_value {
            return Err("Final check failed in Sumcheck protocol");
        }
        println!("Sumcheck 成功");

        Ok(())
    }

    /// Comute initial toatal sum
    fn compute_initial_sum(&self) -> Result<FieldElement, &'static str> {
        let prime = &self.polynomial.terms.values().next().unwrap().prime;
        let mut sum = FieldElement::zero(prime);
        for point in self.iterate_domain(0, vec![])? {
            let value = self.polynomial.multi_poly_evaluate(&point)?;
            sum = sum.add(&value)?;
        }
        Ok(sum)
    }

    /// Compute 
    fn compute_univariate_polynomial(
        &self,
        var_index: usize,
        previous_r_i: &[FieldElement],
    ) -> Result<UnivariatePolynomial, &'static str> {
        let prime = &self.polynomial.terms.values().next().unwrap().prime;
        let mut poly = UnivariatePolynomial::new(vec![FieldElement::zero(prime)]);

        let mut assignment = vec![FieldElement::zero(prime); self.num_variables];

        for (i, r_i) in previous_r_i.iter().enumerate() {
            assignment[i] = r_i.clone();
        }

        let mut vars_to_assign = Vec::new();
        for i in 0..self.num_variables {
            if i != var_index && i >= previous_r_i.len() {
                vars_to_assign.push(i);
            }
        }

        let total_assignments = 1 << vars_to_assign.len();

        for bits in 0..total_assignments {
            for (k, &j) in vars_to_assign.iter().enumerate() {
                let bit = (bits >> k) & 1;
                assignment[j] = FieldElement::from_u32(bit as u32, prime)?;
            }

            let mut term_poly = UnivariatePolynomial::new(vec![FieldElement::zero(prime)]);

            for (exponents, coeff) in &self.polynomial.terms {
                let mut term_coeff = coeff.clone();
                let mut degree = 0u32;

                for (i, &exp) in exponents.iter().enumerate() {
                    if exp == 0 {
                        continue;
                    }
                    if i == var_index {
                        degree = exp;
                    } else {
                        let x_i = &assignment[i];
                        let x_i_pow = x_i.pow(exp.try_into().unwrap())?;
                        term_coeff = term_coeff.mul(&x_i_pow)?;
                    }
                }

                let mut term_poly_coeffs = vec![FieldElement::zero(prime); (degree as usize) + 1];
                term_poly_coeffs[degree as usize] = term_coeff;
                let mono_poly = UnivariatePolynomial::new(term_poly_coeffs);

                term_poly = term_poly.uni_poly_add(&mono_poly)?;
            }

            poly = poly.uni_poly_add(&term_poly)?;
        }

        Ok(poly)
    }

    fn iterate_domain(
        &self,
        var_index: usize,
        current_assignment: Vec<FieldElement>,
    ) -> Result<Vec<Vec<FieldElement>>, &'static str> {
        if var_index >= self.num_variables {
            return Ok(vec![current_assignment]);
        }

        let mut results = Vec::new();
        let prime = &self.polynomial.terms.values().next().unwrap().prime;

        // x_i = 0
        let mut next_assignment_zero = current_assignment.clone();
        next_assignment_zero.push(FieldElement::zero(prime));
        let sub_results_zero = self.iterate_domain(var_index + 1, next_assignment_zero)?;
        results.extend(sub_results_zero);

        // x_i = 1
        let mut next_assignment_one = current_assignment;
        next_assignment_one.push(FieldElement::one(prime));
        let sub_results_one = self.iterate_domain(var_index + 1, next_assignment_one)?;
        results.extend(sub_results_one);

        Ok(results)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::FieldElement;
    use crate::multivariate_polynomial::MultivariatePolynomial;
    use num_bigint::BigUint;
    use std::collections::HashMap;

    #[test]
    fn test_sumcheck_protocol_1() {
        // p = 17
        let prime = BigUint::from(17u32);

        // g(x0, x1, x2) = 2x0^3 + x1 + x0 * x2
        let mut terms = HashMap::new();

        // 2x0^3
        terms.insert(
            vec![3, 0, 0],
            FieldElement::from_u32(2u32, &prime).unwrap(),
        );

        // x1
        terms.insert(vec![0, 1, 0], FieldElement::one(&prime));

        // x0 * x2
        terms.insert(vec![1, 0, 1], FieldElement::one(&prime));

        let polynomial = MultivariatePolynomial::new(terms);

        // r0 = 1, r1 = 0, r2 = 1
        let verifier_challenges = vec![
            FieldElement::from_u32(2u32, &prime).unwrap(),
            FieldElement::from_u32(5u32, &prime).unwrap(),
            FieldElement::from_u32(7u32, &prime).unwrap(),
        ];

        let num_variables = 3;
        let sumcheck = Sumcheck::new(&polynomial, num_variables, verifier_challenges);

        let result = sumcheck.execute();
        assert!(result.is_ok());
    }

    #[test]
    fn test_sumcheck_protocol_2() {
        // p = 17
        let prime = BigUint::from(17u32);

        // g(x0, x1, x2) = 3x0 * x1 + 5x2 + 2
        let mut terms = HashMap::new();

        // 3x0 * x1
        terms.insert(
            vec![1, 1, 0],
            FieldElement::from_u32(3u32, &prime).unwrap(),
        );

        // 5x2
        terms.insert(
            vec![0, 0, 1],
            FieldElement::from_u32(5u32, &prime).unwrap(),
        );

        // 2
        terms.insert(
            vec![0, 0, 0],
            FieldElement::from_u32(2u32, &prime).unwrap(),
        );

        let polynomial = MultivariatePolynomial::new(terms);

        // r1 = 4, r2 = 7, r3 = 5
        let verifier_challenges = vec![
            FieldElement::from_u32(4u32, &prime).unwrap(),
            FieldElement::from_u32(7u32, &prime).unwrap(),
            FieldElement::from_u32(5u32, &prime).unwrap(),
        ];

        let num_variables = 3;
        let sumcheck = Sumcheck::new(&polynomial, num_variables, verifier_challenges);

        let result = sumcheck.execute();
        assert!(result.is_ok());
    }

    #[test]
    fn test_sumcheck_protocol_3() {
        // p = 7
        let prime = BigUint::from(7u32);

        // g(x0, x1) = x0 * x1 + x0 + x1 + 1
        let mut terms = HashMap::new();

        // x0 * x1
        terms.insert(vec![1, 1], FieldElement::one(&prime));

        // x0
        terms.insert(vec![1, 0], FieldElement::one(&prime));

        // x1
        terms.insert(vec![0, 1], FieldElement::one(&prime));

        // 1
        terms.insert(vec![0, 0], FieldElement::one(&prime));

        let polynomial = MultivariatePolynomial::new(terms);

        // r1 = 3, r2 = 2
        let verifier_challenges = vec![
            FieldElement::from_u32(3u32, &prime).unwrap(),
            FieldElement::from_u32(2u32, &prime).unwrap(),
        ];

        let num_variables = 2;
        let sumcheck = Sumcheck::new(&polynomial, num_variables, verifier_challenges);

        let result = sumcheck.execute();
        assert!(result.is_ok());
    }


}
