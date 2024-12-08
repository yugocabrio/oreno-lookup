use crate::field::FieldElement;
use crate::multivariate_polynomial::MultivariatePolynomial;
use crate::univariate_polynomial::UnivariatePolynomial;
use num_bigint::BigUint;

pub struct Prover {
    // Multivariate polynomial that takes in boolean inputs; hypercube
    pub multivariate_polynomial: MultivariatePolynomial,
    pub num_variables: usize,
}

impl Prover {
    pub fn new(multivariate_polynomial: MultivariatePolynomial, num_variables: usize) -> Self {
        Self {
            multivariate_polynomial,
            num_variables,
        }
    }
    
    // Prover computes and fixes an initial sum H in F
    pub fn compute_initial_sum(&self) -> Result<FieldElement, &'static str> {
        let prime = self.multivariate_polynomial.get_prime()?;
        let mut sum = FieldElement::zero(&prime);

        // Generates all possible assignments for the variables from {0,1}^v
        let all_assignments = generate_all_assignments(self.num_variables, &prime)?;
        for assignment in all_assignments {
            // Evaluates the multivariate polynomial at the current assignment and compute toatal initial sum
            let value = self.multivariate_polynomial.multi_poly_evaluate(&assignment)?;
            sum = sum.add(&value)?;
        }
        Ok(sum)
    }

    // In the i-th round, the Prover sends the univariate polynomial and computes partial sums
    /// g_i(X_i) = Σ_{x_{i+1},...,x_n ∈ {0,1}} g(r_1, ..., r_{i-1}, X_i, x_{i+1}, ..., x_n)
    pub fn compute_univariate_polynomial(
        &self,
        var_index: usize, // current variable index
        challenges: &[FieldElement],
    ) -> Result<UnivariatePolynomial, &'static str> {
        let prime = self.multivariate_polynomial.get_prime()?;

        // Initialize the univariate polynomial
        let mut univariate_polynomial =
            UnivariatePolynomial::new(vec![FieldElement::zero(&prime)]);

        // Prepare an assignment vector(zero)
        let mut assignment = vec![FieldElement::zero(&prime); self.num_variables];
        for (i, challenge) in challenges.iter().enumerate() {
            assignment[i] = challenge.clone();
        }

        // Identify the variables to sum over (excluding var_index and fixed challenges)
        let vars_to_sum_over: Vec<usize> = (0..self.num_variables)
            .filter(|&i| i != var_index && i >= challenges.len())
            .collect();

        println!("vars_to_sum_over: {:?}", vars_to_sum_over);

        // Generate all possible assignments that are used for the specific round
        let sum_over_assignments = generate_all_assignments(vars_to_sum_over.len(), &prime)?;
        println!("sum_over_assignments: {:?}", sum_over_assignments);
        for sum_assignment in sum_over_assignments {
            // Update the assignment vector with the current sum assignment
            for (k, &var_idx) in vars_to_sum_over.iter().enumerate() {
                assignment[var_idx] = sum_assignment[k].clone();
            }

            println!("assignment: {:?}", assignment);

            // Compute the univariate polynomial term for the current assignment
            let mut term_polynomial =
                UnivariatePolynomial::new(vec![FieldElement::zero(&prime)]);

            // Iterate over each term in the multivariate polynomial
            for (exponents, coeff) in &self.multivariate_polynomial.terms {
                let mut term_coeff = coeff.clone();
                let mut degree = 0u32;

                // Process each variable in the term
                for (i, &exp) in exponents.iter().enumerate() {
                    if exp == 0 {
                        continue;
                    }
                    if i == var_index {
                        // The current variable remains as X_i in the univariate polynomial
                        degree += exp;
                    } else {
                        // Fixed variables are substituted with their values
                        let x_i = &assignment[i];
                        let x_i_pow = x_i.pow(exp)?;
                        term_coeff = term_coeff.mul(&x_i_pow)?;
                    }
                }

                // Create a monomial for the univariate polynomial
                let mut monomial_coeffs =
                    vec![FieldElement::zero(&prime); (degree as usize) + 1];
                monomial_coeffs[degree as usize] = term_coeff;
                let monomial = UnivariatePolynomial::new(monomial_coeffs);

                // Add the monomial to the term polynomial
                term_polynomial = term_polynomial.uni_poly_add(&monomial)?;
            }

            // Sum the computed univariate polynomial
            univariate_polynomial = univariate_polynomial.uni_poly_add(&term_polynomial)?;
        }

        Ok(univariate_polynomial)
    }

    /// Evaluates the polynomial at the provided challenges g(r_1, ..., r_v).
    pub fn evaluate_multi_variate_polynomial(
        &self,
        challenges: &[FieldElement],
    ) -> Result<FieldElement, &'static str> {
        self.multivariate_polynomial.multi_poly_evaluate(challenges)
    }
}

pub struct Verifier {
    num_variables: usize,
    challenges: Vec<FieldElement>,
    claimed_sums: Vec<FieldElement>,
}

impl Verifier {
    pub fn new(num_variables: usize) -> Self {
        Self {
            num_variables,
            challenges: Vec::new(),
            claimed_sums: Vec::new(),
        }
    }

    // Records a new challenge r_i
    pub fn add_challenge(&mut self, challenge: FieldElement) {
        self.challenges.push(challenge);
    }

    // Records claimed sum g_i(r_i) provided by the Prover
    pub fn add_claimed_sum(&mut self, claimed_sum: FieldElement) {
        self.claimed_sums.push(claimed_sum);
    }

    // Verifier checks that g_{i-1}(r_{i-1}) = g_i(0) + g_i(1)
    pub fn check_univariate_polynomial(
        &self,
        univariate_polynomial: &UnivariatePolynomial,
        previous_claimed_sum: &FieldElement,
    ) -> Result<(), &'static str> {
        let prime = &univariate_polynomial.coefficients[0].prime;

        // Evaluate the univariate polynomial at X_i = 0 and X_i = 1
        let g_0 = univariate_polynomial.uni_poly_evaluate(&FieldElement::zero(prime))?;
        let g_1 = univariate_polynomial.uni_poly_evaluate(&FieldElement::one(prime))?;
        let g_total = g_0.add(&g_1)?;

        // Check that the sum equals the previous claimed sum
        if g_total != *previous_claimed_sum {
            return Err("Verification failed: g(0) + g(1) != previous claimed sum");
        }

        // TODO: degree should be checked?

        Ok(())
    }
}

pub struct Sumcheck {
    pub prover: Prover,
    pub verifier: Verifier,
    pub challenges: Vec<FieldElement>,
    // This is used in the logup protocol
    pub require_initial_sum_zero: bool,
}

impl Sumcheck {
    pub fn new(
        prover: Prover,
        verifier: Verifier,
        challenges: Vec<FieldElement>,
    ) -> Self {
        Self {
            prover,
            verifier,
            challenges,
            require_initial_sum_zero: false,
        }
    }

    pub fn execute_sumcheck_protocol(&mut self) -> Result<(), &'static str> {
        // let prime = self.prover.polynomial.get_prime()?;

        // Prover computes the initial sum H and sends it to the Verifier
        let mut claimed_sum = self.prover.compute_initial_sum()?;
        println!("Initial sum (H): {:?}", claimed_sum.num);

        // This is used in the logup protocol
        if self.require_initial_sum_zero {
            let prime = self.prover.multivariate_polynomial.get_prime()?;
            let zero = FieldElement::zero(&prime);
            if claimed_sum != zero {
                panic!("Initial sum H is not zero! This is required by the logup protocol.");
            }
        }

        // For each round
        for (round, challenge) in self.challenges.iter().enumerate() {
            // Prover computes the univariate polynomial g_i(X_i) and sends it to the Verifier
            let univariate_poly = self
                .prover
                .compute_univariate_polynomial(round, &self.verifier.challenges)?;

            println!(
                "Round {}: Univariate polynomial g{}(X{}) = {:?}",
                round + 1,
                round + 1,
                round + 1,
                univariate_poly
                    .coefficients
                    .iter()
                    .map(|coeff| coeff.num.clone())
                    .collect::<Vec<BigUint>>()
            );

            // Verifier checks the polynomial and the claimed sum
            self.verifier
                .check_univariate_polynomial(&univariate_poly, &claimed_sum)?;

            // Verifier adds the challenge r_i
            self.verifier.add_challenge(challenge.clone());
            println!(
                "Verifier's challenge r{} = {:?}",
                round + 1,
                challenge.num
            );

            // Updates the claimed sum
            claimed_sum = univariate_poly.uni_poly_evaluate(challenge)?;
            self.verifier.add_claimed_sum(claimed_sum.clone());
            println!("Updated claimed sum (H): {:?}", claimed_sum.num);
        }

        // Final round
        // Verifier evaluates multivariate polynomial: g(r_1, ..., r_v).
        let final_evaluation = self
            .prover
            .evaluate_multi_variate_polynomial(&self.verifier.challenges)?;
        println!(
            "Final evaluation of g(r1,...,rn): {:?}",
            final_evaluation.num
        );

        // Verifier checks final evaluation and final claimed sum
        if claimed_sum != final_evaluation {
            return Err("Final verification failed");
        }

        println!("Sum-Check protocol succeeded");
        Ok(())
    }
}

// Generates all possible assignments for variables in {0,1}^n
pub fn generate_all_assignments(
    num_vars: usize,
    prime: &BigUint,
) -> Result<Vec<Vec<FieldElement>>, &'static str> {
    let num_assignments = 1 << num_vars;
    let mut assignments = Vec::with_capacity(num_assignments);

    // Iterate over all integers from 0 to 2^n - 1 to generate binary assignments
    for i in 0..num_assignments {
        let mut assignment = Vec::with_capacity(num_vars);
        for j in (0..num_vars).rev() {
            let bit = (i >> j) & 1;
            assignment.push(FieldElement::from_u32(bit as u32, prime)?);
        }
        assignments.push(assignment);
    }
    Ok(assignments)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::FieldElement;
    use crate::multivariate_polynomial::MultivariatePolynomial;
    use num_bigint::BigUint;
    use std::collections::HashMap;

    #[test]
    fn test_compute_initial_sum() {
        // Prime field P17
        let prime = BigUint::from(17u32);

        // g(x0, x1) = x0 * x1 + x0 + x1 + 1
        let mut terms = HashMap::new();

        // x0 * x1
        terms.insert(vec![1, 1], FieldElement::one(&prime));

        // x0
        terms.insert(vec![1, 0], FieldElement::one(&prime));

        // x1
        terms.insert(vec![0, 1], FieldElement::one(&prime));

        // Constant term 1
        terms.insert(vec![0, 0], FieldElement::one(&prime));

        let polynomial = MultivariatePolynomial::new(terms);

        let num_variables = 2;
        let prover = Prover::new(polynomial, num_variables);

        // g(0,0) + g(0,1) + g(1,0) + g(1,1) = 1 + 2 + 2 + 4 = 9 mod 17
        let initial_sum = prover.compute_initial_sum().unwrap();
        let expected_sum = FieldElement::from_u32(9u32, &prime).unwrap();

        assert_eq!(initial_sum, expected_sum);
    }

    #[test]
    fn test_compute_univariate_polynomial() {
        // Prime field P17
        let prime = BigUint::from(17u32);

        // g(x0, x1) = x0 * x1 + x0 + x1 + 1
        let mut terms = HashMap::new();

        // x0 * x1
        terms.insert(vec![1, 1], FieldElement::one(&prime));

        // x0
        terms.insert(vec![1, 0], FieldElement::one(&prime));

        // x1
        terms.insert(vec![0, 1], FieldElement::one(&prime));

        // Constant term 1
        terms.insert(vec![0, 0], FieldElement::one(&prime));

        let polynomial = MultivariatePolynomial::new(terms);

        let num_variables = 2;
        let prover = Prover::new(polynomial, num_variables);

        // Index 0
        let var_index = 0;
        let previous_challenges = vec![];

        let univariate_poly = prover
            .compute_univariate_polynomial(var_index, &previous_challenges)
            .unwrap();

        // g1(X1) = 3 + 3X1
        let expected_coeffs = vec![
            FieldElement::from_u32(3u32, &prime).unwrap(),
            FieldElement::from_u32(3u32, &prime).unwrap(),
        ];

        assert_eq!(univariate_poly.coefficients, expected_coeffs);
    }

    #[test]
    fn test_generate_all_assignments() {
        // Prime field P17
        let prime = BigUint::from(17u32);

        let num_vars = 2;

        // Expected assignments: [ [0,0], [0,1], [1,0], [1,1] ]
        let expected_assignments = vec![
            vec![
                FieldElement::from_u32(0, &prime).unwrap(),
                FieldElement::from_u32(0, &prime).unwrap(),
            ],
            vec![
                FieldElement::from_u32(0, &prime).unwrap(),
                FieldElement::from_u32(1, &prime).unwrap(),
            ],
            vec![
                FieldElement::from_u32(1, &prime).unwrap(),
                FieldElement::from_u32(0, &prime).unwrap(),
            ],
            vec![
                FieldElement::from_u32(1, &prime).unwrap(),
                FieldElement::from_u32(1, &prime).unwrap(),
            ],
        ];

        let generated_assignments = generate_all_assignments(num_vars, &prime)
            .expect("Failed to generate assignments");

        assert_eq!(generated_assignments, expected_assignments);
    }

    #[test]
    fn test_sumcheck_protocol_1() {
        // Prime field P17
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

        // r0 = 2, r1 = 5, r2 = 7
        let verifier_challenges = vec![
            FieldElement::from_u32(2u32, &prime).unwrap(),
            FieldElement::from_u32(5u32, &prime).unwrap(),
            FieldElement::from_u32(7u32, &prime).unwrap(),
        ];

        let num_variables = 3;
        let prover = Prover::new(polynomial, num_variables);
        let verifier = Verifier::new(num_variables);

        let mut sumcheck = Sumcheck::new(prover, verifier, verifier_challenges);

        let result = sumcheck.execute_sumcheck_protocol();
        assert!(result.is_ok());
    }

    #[test]
    fn test_sumcheck_protocol_2() {
        // Prime field P17
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

        // Constant term 2
        terms.insert(
            vec![0, 0, 0],
            FieldElement::from_u32(2u32, &prime).unwrap(),
        );

        let polynomial = MultivariatePolynomial::new(terms);

        // r0 = 4, r1 = 7, r2 = 5
        let verifier_challenges = vec![
            FieldElement::from_u32(4u32, &prime).unwrap(),
            FieldElement::from_u32(7u32, &prime).unwrap(),
            FieldElement::from_u32(5u32, &prime).unwrap(),
        ];

        let num_variables = 3;
        let prover = Prover::new(polynomial, num_variables);
        let verifier = Verifier::new(num_variables);

        let mut sumcheck = Sumcheck::new(prover, verifier, verifier_challenges);

        let result = sumcheck.execute_sumcheck_protocol();
        assert!(result.is_ok());
    }

    #[test]
    fn test_sumcheck_protocol_3() {
        // Prime field P7
        let prime = BigUint::from(7u32);

        // g(x0, x1) = x0 * x1 + x0 + x1 + 1
        let mut terms = HashMap::new();

        // x0 * x1
        terms.insert(vec![1, 1], FieldElement::one(&prime));

        // x0
        terms.insert(vec![1, 0], FieldElement::one(&prime));

        // x1
        terms.insert(vec![0, 1], FieldElement::one(&prime));

        // Constant term 1
        terms.insert(vec![0, 0], FieldElement::one(&prime));

        let polynomial = MultivariatePolynomial::new(terms);

        // r0 = 3, r1 = 2
        let verifier_challenges = vec![
            FieldElement::from_u32(3u32, &prime).unwrap(),
            FieldElement::from_u32(2u32, &prime).unwrap(),
        ];

        let num_variables = 2;
        let prover = Prover::new(polynomial, num_variables);
        let verifier = Verifier::new(num_variables);

        let mut sumcheck = Sumcheck::new(prover, verifier, verifier_challenges);

        let result = sumcheck.execute_sumcheck_protocol();
        assert!(result.is_ok());
    }
}
