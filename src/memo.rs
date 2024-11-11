use crate::field::FieldElement;
use crate::multivariate_polynomial::MultivariatePolynomial;
use num_bigint::BigUint;
use std::collections::HashMap;

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
            // Prover computes univariate polynomial coefficients S_0 and S_1
            let (s_0, s_1) = self.compute_univariate_coefficients(i, &previous_r_i)?;
            println!(
                "ラウンド {}: 一変数多項式係数 g{}(X{}) = S_0: {:?}, S_1: {:?}",
                i + 1,
                i + 1,
                i + 1,
                s_0,
                s_1
            );

            // Verifier checks that s_0 + s_1 == claimed_sum
            let s_total = s_0.add(&s_1)?;
            if claimed_sum != s_total {
                return Err("Consistency check failed in Sumcheck protocol");
            }
            println!(
                "検証者は一致性を確認しました。S_0 + S_1 = claimed_sum = {:?}",
                claimed_sum
            );

            // Verifier provides challenge (from pre-defined list)
            let r_i = self
                .verifier_challenges
                .get(i)
                .ok_or("Not enough challenges provided")?;
            println!("検証者が選んだランダムチャレンジ r{} = {:?}", i + 1, r_i);

            // Prover and Verifier update claimed sum
            claimed_sum = s_0.mul(&FieldElement::one(prime).sub(r_i)?)?.add(&s_1.mul(r_i)?)?;
            println!("更新された総和：{:?}", claimed_sum);

            previous_r_i.push(r_i.clone());
        }

        // Final check
        let final_value = self.polynomial.multi_poly_evaluate(&previous_r_i)?;
        println!("最終評価点での多項式値：{:?}", final_value);
        if claimed_sum != final_value {
            return Err("Final check failed in Sumcheck protocol");
        }
        println!("Sumcheck プロトコルが成功しました。");

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

    /// 一変数多項式の係数 S_0, S_1 を計算
    fn compute_univariate_coefficients(
        &self,
        var_index: usize,
        previous_r_i: &[FieldElement],
    ) -> Result<(FieldElement, FieldElement), &'static str> {
        // x_i = 0 の場合の部分総和 S_0
        let mut assignment_zero = previous_r_i.to_vec();
        assignment_zero.push(FieldElement::zero(
            &self.polynomial.terms.values().next().unwrap().prime,
        ));
        let s_0 = self.compute_partial_sum(var_index + 1, &assignment_zero)?;

        // x_i = 1 の場合の部分総和 S_1
        let mut assignment_one = previous_r_i.to_vec();
        assignment_one.push(FieldElement::one(
            &self.polynomial.terms.values().next().unwrap().prime,
        ));
        let s_1 = self.compute_partial_sum(var_index + 1, &assignment_one)?;

        Ok((s_0, s_1))
    }

    /// 部分総和を計算
    fn compute_partial_sum(
        &self,
        var_start_index: usize,
        assignment: &[FieldElement],
    ) -> Result<FieldElement, &'static str> {
        if var_start_index >= self.num_variables {
            return self.polynomial.multi_poly_evaluate(assignment);
        }

        let prime = &self.polynomial.terms.values().next().unwrap().prime;
        let mut sum = FieldElement::zero(prime);

        // x_i = 0 の場合
        let mut assignment_zero = assignment.to_vec();
        assignment_zero.push(FieldElement::zero(prime));
        let partial_sum_zero =
            self.compute_partial_sum(var_start_index + 1, &assignment_zero)?;
        sum = sum.add(&partial_sum_zero)?;

        // x_i = 1 の場合
        let mut assignment_one = assignment.to_vec();
        assignment_one.push(FieldElement::one(prime));
        let partial_sum_one = self.compute_partial_sum(var_start_index + 1, &assignment_one)?;
        sum = sum.add(&partial_sum_one)?;

        Ok(sum)
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
            FieldElement::from_u32(1u32, &prime).unwrap(), // r0 = 1
            FieldElement::from_u32(0u32, &prime).unwrap(), // r1 = 0
            FieldElement::from_u32(1u32, &prime).unwrap(), // r2 = 1
        ];

        let num_variables = 3;
        let sumcheck = Sumcheck::new(&polynomial, num_variables, verifier_challenges);

        let result = sumcheck.execute();
        assert!(result.is_ok());
    }
}
