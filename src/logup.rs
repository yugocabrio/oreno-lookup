use std::collections::HashMap;
use crate::field::FieldElement;
use crate::multivariate_polynomial::MultivariatePolynomial;
use crate::sumcheck::{Prover as SumcheckProver, Verifier as SumcheckVerifier, Sumcheck};
use num_bigint::BigUint;

pub struct LogUpProver {
    // Witness columns: w1(X), ..., wM(X)
    witness_columns: Vec<Vec<FieldElement>>,
    // Table values: t(X)
    table_values: Vec<FieldElement>,
    // Multiplicities: m(X)
    multiplicities: Vec<FieldElement>,
    // Prime p
    prime: BigUint,
}

impl LogUpProver {
    pub fn new(
        witness_columns: Vec<Vec<FieldElement>>,
        table_values: Vec<FieldElement>,
        multiplicities: Vec<FieldElement>,
        prime: BigUint,
    ) -> Self {
        Self {
            witness_columns,
            table_values,
            multiplicities,
            prime,
        }
    }

    /// Constructs the multivariate polynomial P(x) using Lagrange basis.
    pub fn construct_polynomial(
        &self,
        alpha: &FieldElement,
    ) -> Result<MultivariatePolynomial, &'static str> {
        // Number of rows (assignments)
        let num_rows = self.witness_columns[0].len();
        // Number of variables n such that 2^n = num_rows
        let n = (num_rows as f64).log2() as usize;

        let all_assignments = generate_all_assignments(n, &self.prime)?;

        // Verify the number of assignments matches
        if all_assignments.len() != num_rows {
            return Err("Number of assignments does not match number of rows");
        }

        // Initialize P(x) as the zero polynomial
        let mut P = MultivariatePolynomial::new(HashMap::new());

        println!("Constructing polynomial P(x) for all x in H_n");
        println!("Number of variables (n): {}", n);
        println!("Number of assignments: {}", all_assignments.len());

        // Iterate over each assignment x ∈ H_n
        for (idx, x_assignment) in all_assignments.iter().enumerate() {
            // Retrieve m(x) and t(x)
            let m_x = self.multiplicities[idx].clone();
            let t_x = self.table_values[idx].clone();

            // Retrieve witness values w_i(x)
            let w_values: Vec<FieldElement> = self.witness_columns.iter().map(|w_col| w_col[idx].clone()).collect();

            // Compute f(x) = m(x)/(alpha - t(x)) - sum_i 1/(alpha - w_i(x))
            let alpha_sub_t_x = alpha.sub(&t_x)?;
            let numerator = m_x.div(&alpha_sub_t_x)?;

            let mut denominator_sum = FieldElement::zero(&self.prime);
            for w_i_x in &w_values {
                let alpha_sub_w_i_x = alpha.sub(w_i_x)?;
                let term = FieldElement::one(&self.prime).div(&alpha_sub_w_i_x)?;
                denominator_sum = denominator_sum.add(&term)?;
            }

            let f_x = numerator.sub(&denominator_sum)?;

            // Only include non-zero terms
            if !f_x.is_zero() {
                // Construct Lagrange basis polynomial L_x(x)
                let l_x = construct_lagrange_basis(x_assignment, n)?;

                // Multiply L_x(x) by f(x): P(x) += f(x) * L_x(x)
                let f_x_l_x = l_x.multi_poly_scale(&f_x)?;
                P = P.multi_poly_add(&f_x_l_x)?;
                
                // Debug information
                println!("x_assignment[{}]: {:?}", idx, x_assignment.iter().map(|e| e.num.clone()).collect::<Vec<BigUint>>());
                println!("f(x): {}", f_x.num);
                println!("L_x(x): {:?}", l_x);
                println!("f(x) * L_x(x): {:?}", f_x_l_x);
                println!("-----------------------------");
            }
        }

        println!("polynomial P(x): {:?}", P);
        Ok(P)
    }

    /// Executes the Sumcheck protocol.
    pub fn execute_sumcheck(&self, alpha: &FieldElement) -> Result<(), &'static str> {
        // Construct the multivariate polynomial P(x)
        let P = self.construct_polynomial(alpha)?;

        // Number of variables
        let num_rows = self.witness_columns[0].len();
        let n = (num_rows as f64).log2() as usize;

        // Initialize Sumcheck Prover and Verifier
        let sumcheck_prover = SumcheckProver::new(P.clone(), n);
        let sumcheck_verifier = SumcheckVerifier::new(n);


        // Define challenges (here using fixed challenges; in practice, these should be random)
        let challenges = (0..n)
            .map(|i| FieldElement::from_u32((i % 2) as u32, &self.prime))
            .collect::<Result<Vec<FieldElement>, _>>()?;

        println!("Challenges (Verifier's random values): {:?}", challenges.iter().map(|e| e.num.clone()).collect::<Vec<BigUint>>());

        // Initialize the Sumcheck protocol
        let mut sumcheck = Sumcheck::new(sumcheck_prover, sumcheck_verifier, challenges);

        // In the logup protocol, initial sum H must be 0.
        sumcheck.require_initial_sum_zero = true;

        // Execute the Sumcheck protocol
        sumcheck.execute_sumcheck_protocol()
    }
}

// Constructs the Lagrange basis polynomial L_x(x) for a given assignment.
// L_x(x) = Π_{i=1}^n [ (x_i)^a_i * (1 - x_i)^(1 - a_i) ]
fn construct_lagrange_basis(x_assignment: &[FieldElement], n: usize) -> Result<MultivariatePolynomial, &'static str> {
    let prime = &x_assignment[0].prime;
    let mut basis = MultivariatePolynomial::new(HashMap::new());

    // Initialize basis as the constant polynomial 1
    basis.terms.insert(vec![0; n], FieldElement::one(prime).clone());

    for (i, bit) in x_assignment.iter().enumerate() {
        let mut term = MultivariatePolynomial::new(HashMap::new());

        // Compare BigUint with BigUint::from(1u32)
        if bit.num == BigUint::from(1u32) {
            // Multiply by x_i
            let mut exponents = vec![0; n];
            exponents[i] = 1;
            term.terms.insert(exponents, FieldElement::one(prime).clone());
        } else {
            // Multiply by (1 - x_i) = 1 + (-1) * x_i
            let mut exponents_zero = vec![0; n];
            term.terms.insert(exponents_zero.clone(), FieldElement::one(prime).clone());

            let mut exponents_one = vec![0; n];
            exponents_one[i] = 1;
            // Use negate() to get -1
            let neg_one = FieldElement::one(prime).negate()?;
            term.terms.insert(exponents_one.clone(), neg_one);
        }

        // Multiply the current basis by the term
        basis = basis.multi_poly_mul(&term)?;
    }

    Ok(basis)
}

// Generates all possible assignments for variables in {0,1}^n
pub fn generate_all_assignments(
    num_vars: usize,
    prime: &BigUint,
) -> Result<Vec<Vec<FieldElement>>, &'static str> {
    let num_assignments = 1 << num_vars;
    let mut assignments = Vec::with_capacity(num_assignments);

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
    use num_bigint::BigUint;

    #[test]
    fn test_logup_protocol_test_1() {
        let prime = BigUint::from(97u32);
        let w1 = vec![
            FieldElement::from_u32(1, &prime).unwrap(),
            FieldElement::from_u32(2, &prime).unwrap(),
            FieldElement::from_u32(3, &prime).unwrap(),
            FieldElement::from_u32(1, &prime).unwrap(),
        ];
        let w2 = vec![
            FieldElement::from_u32(2, &prime).unwrap(),
            FieldElement::from_u32(1, &prime).unwrap(),
            FieldElement::from_u32(4, &prime).unwrap(),
            FieldElement::from_u32(2, &prime).unwrap(),
        ];

        let t = vec![
            FieldElement::from_u32(1, &prime).unwrap(),
            FieldElement::from_u32(2, &prime).unwrap(),
            FieldElement::from_u32(3, &prime).unwrap(),
            FieldElement::from_u32(4, &prime).unwrap(),
        ];

        let m = vec![
            FieldElement::from_u32(3, &prime).unwrap(),
            FieldElement::from_u32(3, &prime).unwrap(),
            FieldElement::from_u32(1, &prime).unwrap(),
            FieldElement::from_u32(1, &prime).unwrap(),
        ];

        let prover = LogUpProver::new(vec![w1.clone(), w2.clone()], t.clone(), m.clone(), prime.clone());

        let alpha = FieldElement::from_u32(5, &prime).unwrap();

        let result = prover.execute_sumcheck(&alpha);
        assert!(result.is_ok());
    }
}
