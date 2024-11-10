use crate::field::FieldElement;
use num_bigint::BigUint;
use num_traits::FromPrimitive;

pub fn compute_F(
    f_values: &[u32],
    t_values: &[u32],
    beta: &FieldElement,
    gamma: &FieldElement,
    prime: &BigUint,
) -> Result<FieldElement, &'static str> {
    let one = FieldElement::one(prime);
    let one_plus_beta = one.add(beta)?;
    let n = f_values.len();

    // (1 + β)^n
    let mut result = one_plus_beta.pow(n as u32)?;

    // Π {i∈[n]} (γ + f_i)
    for &fi in f_values {
        let fi_fe = FieldElement::new(BigUint::from(fi), prime.clone())?;
        result = result.mul(&gamma.add(&fi_fe)?)?;
    }

    // ∏_{i∈[d−1]} (γ(1 + β) + t_i + β t_{i+1})
    let gamma_one_plus_beta = gamma.mul(&one_plus_beta)?;
    for i in 0..(t_values.len() - 1) {
        let t_i = FieldElement::new(BigUint::from(t_values[i]), prime.clone())?;
        let t_i_plus1 = FieldElement::new(BigUint::from(t_values[i + 1]), prime.clone())?;
        let term = gamma_one_plus_beta.add(&t_i)?.add(&t_i_plus1.mul(beta)?)?;
        result = result.mul(&term)?;
    }

    Ok(result)
}

pub fn compute_G(
    s_values: &[u32],
    beta: &FieldElement,
    gamma: &FieldElement,
    prime: &BigUint,
) -> Result<FieldElement, &'static str> {
    let one = FieldElement::one(prime);
    let one_plus_beta = one.add(beta)?;
    let gamma_one_plus_beta = gamma.mul(&one_plus_beta)?;
    let mut result = FieldElement::one(prime);

    // ∏_{i∈[2n]} (γ(1 + β) + s_i + β s_{i+1})
    for i in 0..(s_values.len() - 1) {
        let s_i = FieldElement::new(BigUint::from(s_values[i]), prime.clone())?;
        let s_i_plus1 = FieldElement::new(BigUint::from(s_values[i + 1]), prime.clone())?;
        let term = gamma_one_plus_beta.add(&s_i)?.add(&s_i_plus1.mul(beta)?)?;
        result = result.mul(&term)?;
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::FieldElement;
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;

    #[test]
    fn test_naive_plookup() {
        let prime = BigUint::from_u32(31).unwrap();

        let f_values = vec![1, 1, 4, 8];
        let t_values = vec![1, 4, 8];

        // s_values
        let mut s_values = f_values.clone();
        s_values.extend(t_values.clone());
        s_values.sort();

        let beta = FieldElement::new(BigUint::from_u32(3).unwrap(), prime.clone()).unwrap();
        let gamma = FieldElement::new(BigUint::from_u32(5).unwrap(), prime.clone()).unwrap();

        let F_result = compute_F(&f_values, &t_values, &beta, &gamma, &prime).unwrap();
        let G_result = compute_G(&s_values, &beta, &gamma, &prime).unwrap();

        println!("F(β, γ) = {}", F_result.num);
        println!("G(β, γ) = {}", G_result.num);

        assert_eq!(F_result, G_result);
    }
}
