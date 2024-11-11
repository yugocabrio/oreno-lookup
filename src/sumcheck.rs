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
        println!("初期総和（H）：{:?}", claimed_sum);
        let mut previous_r_i = Vec::new();

        for i in 0..self.num_variables {
            // Prover computes univariate polynomial
            let univariate_poly = self.compute_univariate_polynomial(i, &previous_r_i)?;
            println!(
                "ラウンド {}: 一変数多項式 g{}(x{}) = {:?}",
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
                "一変数多項式の評価: S_0 = {:?}, S_1 = {:?}",
                s_0, s_1
            );

            // Verifier checks that s_0 + s_1 == claimed_sum
            let s_total = s_0.add(&s_1)?;
            if claimed_sum != s_total {
                return Err("Consistency check failed in Sumcheck protocol");
            }
            println!(
                "S_0 + S_1 = claimed_sum = {:?}",
                claimed_sum
            );

            // Verifier provides challenge (from pre-defined list)
            let r_i = self
                .verifier_challenges
                .get(i)
                .ok_or("Not enough challenges provided")?;
            println!("ランダムチャレンジ r{} = {:?}", i + 1, r_i);

            // Prover and Verifier update claimed sum
            claimed_sum = univariate_poly.uni_poly_evaluate(r_i)?;
            println!("更新された総和：{:?}", claimed_sum);

            previous_r_i.push(r_i.clone());
        }

        // Final check
        let final_value = self.polynomial.multi_poly_evaluate(&previous_r_i)?;
        println!("最終評価点での多項式値：{:?}", final_value);
        if claimed_sum != final_value {
            return Err("Final check failed in Sumcheck protocol");
        }
        println!("Sumcheck 成功");

        Ok(())
    }

    /// 初期の総和を計算
    fn compute_initial_sum(&self) -> Result<FieldElement, &'static str> {
        let prime = &self.polynomial.terms.values().next().unwrap().prime;
        let mut sum = FieldElement::zero(prime);
        for point in self.iterate_domain(0, vec![])? {
            let value = self.polynomial.multi_poly_evaluate(&point)?;
            sum = sum.add(&value)?;
        }
        Ok(sum)
    }

    /// 一変数多項式をシンボリックに計算
    fn compute_univariate_polynomial(
        &self,
        var_index: usize,
        previous_r_i: &[FieldElement],
    ) -> Result<UnivariatePolynomial, &'static str> {
        let prime = &self.polynomial.terms.values().next().unwrap().prime;
        let mut poly = UnivariatePolynomial::new(vec![FieldElement::zero(prime)]);

        // assignment を self.num_variables の長さで初期化
        let mut assignment = vec![FieldElement::zero(prime); self.num_variables];

        // previous_r_i の値を設定
        for (i, r_i) in previous_r_i.iter().enumerate() {
            assignment[i] = r_i.clone();
        }

        // xi の位置はそのまま（変数として残す）
        // その他の変数のインデックスを取得
        let mut vars_to_assign = Vec::new();
        for i in 0..self.num_variables {
            if i != var_index && i >= previous_r_i.len() {
                vars_to_assign.push(i);
            }
        }

        let total_assignments = 1 << vars_to_assign.len();

        for bits in 0..total_assignments {
            // 他の変数を割り当て
            for (k, &j) in vars_to_assign.iter().enumerate() {
                let bit = (bits >> k) & 1;
                assignment[j] = FieldElement::from_u32(bit as u32, prime)?;
            }

            // 各項を計算
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

                // 多項式に項を追加
                let mut term_poly_coeffs = vec![FieldElement::zero(prime); (degree as usize) + 1];
                term_poly_coeffs[degree as usize] = term_coeff;
                let mono_poly = UnivariatePolynomial::new(term_poly_coeffs);

                term_poly = term_poly.uni_poly_add(&mono_poly)?;
            }

            poly = poly.uni_poly_add(&term_poly)?;
        }

        Ok(poly)
    }

    /// 領域を反復処理するためのヘルパー関数
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

        // x_i = 0 の場合
        let mut next_assignment_zero = current_assignment.clone();
        next_assignment_zero.push(FieldElement::zero(prime));
        let sub_results_zero = self.iterate_domain(var_index + 1, next_assignment_zero)?;
        results.extend(sub_results_zero);

        // x_i = 1 の場合
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
    fn test_sumcheck_protocol() {
        // フィールドの素数を設定 (p = 17)
        let prime = BigUint::from(17u32);

        // 多変数多項式を定義：g(x0, x1, x2) = 2x0^3 + x1 + x0 * x2
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

        // 検証者のチャレンジを設定
        // r0 = 1, r1 = 0, r2 = 1
        let verifier_challenges = vec![
            FieldElement::from_u32(2u32, &prime).unwrap(), // r0 = 1
            FieldElement::from_u32(5u32, &prime).unwrap(), // r1 = 0
            FieldElement::from_u32(7u32, &prime).unwrap(), // r2 = 1
        ];

        let num_variables = 3;
        let sumcheck = Sumcheck::new(&polynomial, num_variables, verifier_challenges);

        let result = sumcheck.execute();
        assert!(result.is_ok());
    }
}
