use num_bigint::BigUint;
use num_traits::{FromPrimitive, ToPrimitive, Zero, One};
use crate::field::FieldElement;
use crate::univariate_polynomial::UnivariatePolynomial;

pub fn plookup_protocol(
    // Lookuped elements.
    f_values: Vec<u32>,
    // Table elements.
    t_values: Vec<u32>,
    prime: &BigUint,
) -> Result<bool, &'static str> {
    // Protocol 1:
    // Let s ∈ F^2n+1 be the vector that is (f, t) sorted by t.
    // We represent s by h1, h2 ∈ F<n+1[X] as follows. h_1(g^i) = si for i ∈ [n + 1]; and h_2(g^i) = sn+i for each i ∈ [n + 1].

    // n is the number of elements that is lookuped
    // h1 and h2 is defined as the polynomial whose degree is n+1
    // padding t to be n+1

    // s will be created in the next section
    let (n, n_plus_one, t_padded) = setup_parameters(&f_values, &t_values)?;


    // Protocol 2:
    // P computes the polynomials h1, h2 and sends them to the ideal party I.

    // Prover constructs s , and then polynomial h1,h2
    // Create a multiplicative cyclic subgroup of the required length for lagrange interpolation
    let (p, g, H) = define_subgroup(n_plus_one, prime)?;
    let (s, h1_poly, h2_poly) = construct_s_h1_h2(&f_values, &t_values, &H, &p, n)?;

    // Protocol 3:
    // V chooses random β, γ ∈ F and sends them to P.
    let beta = FieldElement::new(BigUint::from_u32(7u32).unwrap(), p.clone()).unwrap();
    let gamma = FieldElement::new(BigUint::from_u32(11u32).unwrap(), p.clone()).unwrap();
    println!("beta = {}, gamma = {}", beta.num, gamma.num);

    // Protocol 4:
    // P computes a polynomial Z ∈ F<n+1[X] that aggregates the value F(β, γ)/G(β, γ)
    let Z_values = compute_Z_values(&f_values, &t_padded, &s, &beta, &gamma, &p, n)?;
    let Z_poly = UnivariatePolynomial::uni_poly_lagrange_interpolation(&H, &Z_values)?;
    println!("Z_poly = {:?}", Z_poly.coefficients);

    // Protocol 5:
    // 5. P sends Z to I.

    // 6. V checks that Z is indeed of the form described above, and that Z(g^n+1) = 1.
    // More precisely, V checks the following identities for all x ∈ H.
    let all_checks_pass = verify_identities(
        &H, &g, n, &f_values, &t_padded, &s, &beta, &gamma,
        &h1_poly, &h2_poly, &Z_poly, &p,
    )?;

    Ok(all_checks_pass)
}

fn setup_parameters(
    f_values: &[u32],
    t_values: &[u32],
) -> Result<(usize, usize, Vec<u32>), &'static str> {
    // n is the number of elements that is lookuped; len of f
    let n = f_values.len();
    // h1 and h2 is defined as the polynomial whose degree is n+1
    let n_plus_one = n + 1;
    println!("n_plus_one: {:?}", n_plus_one);
    let mut t_padded = t_values.to_vec();

    // Padding t to be n+1
    if t_padded.len() <= n {
        let padding = n - t_padded.len() + 1;
        let last_value = *t_padded.last().unwrap();
        t_padded.extend(vec![last_value; padding]);
    }
    println!("t_padded: {:?}", t_padded);

    Ok((n, n_plus_one, t_padded))
}

fn define_subgroup(
    n_plus_one: usize,
    prime: &BigUint,
) -> Result<(BigUint, FieldElement, Vec<FieldElement>), &'static str> {
    let p = prime.clone();

    // Find n+1 subgroups to get n+1 length, the domain for constructing polynomials in Fp, h1, and h2
    let primitive_root = FieldElement::find_primitive_root(n_plus_one, &p)?;
    let g = FieldElement::new(primitive_root.clone(), p.clone())?;

    // Check g^(n+1) is 1
    let g_to_n_plus_one = g.pow(n_plus_one as u32)?;
    assert_eq!(g_to_n_plus_one.num, BigUint::one(), "g^(n+1) must be 1");

    // Multiplicative cyclic subgroup, n+1 elements
    let mut H = Vec::new();
    for i in 1..=n_plus_one {
        let g_i = g.pow(i as u32)?;
        H.push(g_i);
    }
    println!("multiplicative cyclic subgroup H: {:?}", H);

    Ok((p, g, H))
}

fn construct_s_h1_h2(
    f_values: &[u32],
    t_values: &[u32],
    H: &[FieldElement],
    p: &BigUint,
    n: usize,
) -> Result<(Vec<u32>, UnivariatePolynomial, UnivariatePolynomial), &'static str> {
    // Get f and t together and sort to make s
    let mut s = f_values.to_vec();
    s.extend(t_values.to_vec());
    s.sort();
    println!("s: {:?}", s);

    // Padding to make s length 2n+1
    while s.len() < 2 * n + 1 {
        s.push(*s.last().unwrap());
    }
    println!("s after padding: {:?}", s);
    assert_eq!(s.len(), 2 * n + 1);

    // gets elements h1 and h2
    let mut h1_values = Vec::new();
    let mut h2_values = Vec::new();
    for i in 0..=n {
        h1_values.push(s[i]);
        h2_values.push(s[n + i]);
    }

    println!("h1_values: {:?}", h1_values);
    println!("h2_values: {:?}", h2_values);

    let h1_field_values: Vec<FieldElement> = h1_values
        .iter()
        .map(|&val| FieldElement::new(BigUint::from(val), p.clone()).unwrap())
        .collect();

    let h2_field_values: Vec<FieldElement> = h2_values
        .iter()
        .map(|&val| FieldElement::new(BigUint::from(val), p.clone()).unwrap())
        .collect();

    // Constructs h1 and h2 by Lagrange interpolation
    let h1_poly = UnivariatePolynomial::uni_poly_lagrange_interpolation(H, &h1_field_values)?;
    let h2_poly = UnivariatePolynomial::uni_poly_lagrange_interpolation(H, &h2_field_values)?;

    println!("h1_poly: {:?}", h1_poly.coefficients);
    println!("h2_poly: {:?}", h2_poly.coefficients);

    Ok((s, h1_poly, h2_poly))
}

fn compute_Z_values(
    f_values: &[u32],
    t: &[u32],
    s: &[u32],
    beta: &FieldElement,
    gamma: &FieldElement,
    p: &BigUint,
    n: usize,
) -> Result<Vec<FieldElement>, &'static str> {
    let one = FieldElement::one(&p);
    let one_plus_beta = one.add(beta)?;
    let mut Z_values = Vec::new();

    // (a) Z(g) = 1,
    Z_values.push(one.clone());
    println!("Z_values: {:?}", Z_values);

    // (b) For 2 ≤ i ≤ n
    for i in 2..=n {
        // (1 + β) ^ i - 1
        let exponent = (i - 1) as u32;
        let one_plus_beta_pow = one_plus_beta.pow(exponent)?;
        let mut numerator = one_plus_beta_pow;

        // (1 + β) ^ i - 1 * ∏_{j < i} (γ + f_j)
        for j in 0..(i - 1) {
            let gamma_plus_fj = gamma.add(&FieldElement::new(BigUint::from(f_values[j]), p.clone())?)?;
            numerator = numerator.mul(&gamma_plus_fj)?;
        }

        // (1 + β) ^ i - 1 * ∏_{j < i} (γ + f_j) * ∏_{1 ≤ j < i} (γ(1 + β) + t_j + β t_{j+1})
        for j in 1..i {
            let gamma_one_plus_beta = gamma.mul(&one_plus_beta)?;
            let t_j = FieldElement::new(BigUint::from(t[j - 1]), p.clone())?;
            let t_j_plus1 = FieldElement::new(BigUint::from(t[j]), p.clone())?;
            let term = gamma_one_plus_beta.add(&t_j)?.add(&t_j_plus1.mul(beta)?)?;
            numerator = numerator.mul(&term)?;
        }

        let mut denominator = one.clone();

        // ∏_{1 ≤ j < i} (γ(1 + β) + s_j + β s_{j+1})
        for j in 1..i {
            let gamma_one_plus_beta = gamma.mul(&one_plus_beta)?;
            let s_j = FieldElement::new(BigUint::from(s[j - 1]), p.clone())?;
            let s_j_plus1 = FieldElement::new(BigUint::from(s[j]), p.clone())?;
            let term = gamma_one_plus_beta.add(&s_j)?.add(&s_j_plus1.mul(beta)?)?;
            denominator = denominator.mul(&term)?;
        }

        // ∏_{1 ≤ j < i} (γ(1 + β) + s_j + β s_{j+1}) * (γ(1 + β) + s_{n + j} + β s_{n + j + 1})
        for j in 1..i {
            let gamma_one_plus_beta = gamma.mul(&one_plus_beta)?;
            let s_nj = FieldElement::new(BigUint::from(s[n + j - 1]), p.clone())?;
            let s_nj_plus1 = FieldElement::new(BigUint::from(s[n + j]), p.clone())?;
            let term = gamma_one_plus_beta.add(&s_nj)?.add(&s_nj_plus1.mul(beta)?)?;
            denominator = denominator.mul(&term)?;
        }

        let denominator_inv = denominator.inv()?;
        let Z_i = numerator.mul(&denominator_inv)?;
        Z_values.push(Z_i);
    }
    println!("Z_values: {:?}", Z_values);

    // (c) Z(g^n+1) = 1.
    Z_values.push(one.clone());
    println!("Z_values: {:?}", Z_values);

    println!("Z_values length: {}", Z_values.len());

    Ok(Z_values)
}

fn verify_identities(
    H: &[FieldElement],
    g: &FieldElement,
    n: usize,
    f_values: &[u32],
    t: &[u32],
    s: &[u32],
    beta: &FieldElement,
    gamma: &FieldElement,
    h1_poly: &UnivariatePolynomial,
    h2_poly: &UnivariatePolynomial,
    Z_poly: &UnivariatePolynomial,
    p: &BigUint,
) -> Result<bool, &'static str> {
    let one = FieldElement::one(p);
    let one_plus_beta = one.add(beta)?;

    let n_plus_one = n + 1;
    let n_field = FieldElement::new(BigUint::from(n_plus_one as u32), p.clone())?;
    // x = g^1 で 1
    let L1 = UnivariatePolynomial::uni_poly_lagrange_basis(H, 0, &n_field)?;
    // println!("L1: {:?}", L1);
    // x = g^{n+1} で 1
    let Ln1 = UnivariatePolynomial::uni_poly_lagrange_basis(H, n, &n_field)?;
    // println!("L1: {:?}", Ln1);

    let mut all_checks_pass = true;

    // V checks the following identities for all x ∈ H
    for i in 0..H.len() {
        let x = &H[i];
        let x_next = if i + 1 < H.len() {
            &H[i + 1]
        } else {
            &g.pow((n_plus_one + 1) as u32)? // x = g^{n+2}
        };

        // f(x)
        let f_x = if i < f_values.len() {
            FieldElement::new(BigUint::from(f_values[i]), p.clone())?
        } else {
            FieldElement::zero(p)
        };

        // t(x)
        let t_x = if i < t.len() {
            FieldElement::new(BigUint::from(t[i]), p.clone())?
        } else {
            FieldElement::zero(p)
        };

        // t(g · x)
        let t_gx = if i + 1 < t.len() {
            FieldElement::new(BigUint::from(t[i + 1]), p.clone())?
        } else {
            FieldElement::zero(p)
        };

        // h1(x), h1(g · x)
        let h1_x = h1_poly.uni_poly_evaluate(x)?;
        let h1_gx = h1_poly.uni_poly_evaluate(x_next)?;

        // h2(x), h2(g · x)
        let h2_x = h2_poly.uni_poly_evaluate(x)?;
        let h2_gx = h2_poly.uni_poly_evaluate(x_next)?;

        // Z(x), Z(g · x)
        let Z_x = Z_poly.uni_poly_evaluate(x)?;
        let Z_gx = Z_poly.uni_poly_evaluate(x_next)?;

        // L1(x), Ln+1(x)
        let L1_x = L1.uni_poly_evaluate(x)?;
        let Ln1_x = Ln1.uni_poly_evaluate(x)?;

        // (a): L1(x)(Z(x) − 1) = 0.
        if !L1_x.mul(&Z_x.sub(&one)?)?.is_zero() {
            println!("(a) is failed at x = {}", x.num);
            all_checks_pass = false;
        }

        // (d): Ln+1(x)(Z(x) − 1) = 0.
        if !Ln1_x.mul(&Z_x.sub(&one)?)?.is_zero() {
            println!("(d) is failed at x = {}", x.num);
            all_checks_pass = false;
        }

        // (c): Ln+1(x)(h1(x) − h2(g · x)) = 0.
        if !Ln1_x.mul(&h1_x.sub(&h2_gx)?)?.is_zero() {
            println!("(d) is failed at x = {}", x.num);
            all_checks_pass = false;
        }

        // (b):
        // (x − g^n+1)Z(x)(1 + β) · (γ + f(x))(γ(1 + β) + t(x) + βt(g · x))
        let x_minus_g_n1 = x.sub(&g.pow(n_plus_one as u32)?)?;
        let lhs_b = x_minus_g_n1
            .mul(&Z_x)?
            .mul(&one_plus_beta)?
            .mul(&gamma.add(&f_x)?)?
            .mul(&gamma.mul(&one_plus_beta)?.add(&t_x)?.add(&t_gx.mul(beta)?)?)?;

        //  (x − g^n+1)Z(g·x) (γ(1 + β) + h1(x) + βh1(g·x)) (γ(1 + β)+h2(x) + βh2(g·x))
        let rhs_b = x_minus_g_n1
            .mul(&Z_gx)?
            .mul(&gamma.mul(&one_plus_beta)?.add(&h1_x)?.add(&h1_gx.mul(beta)?)?)?
            .mul(&gamma.mul(&one_plus_beta)?.add(&h2_x)?.add(&h2_gx.mul(beta)?)?)?;

        if lhs_b != rhs_b {
            println!("(b) is failed at x = {}", x.num);
            println!("lhs_b: {}", lhs_b.num);
            println!("rhs_b: {}", rhs_b.num);
            all_checks_pass = false;
        }
    }

    Ok(all_checks_pass)
}


#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    #[test]
    fn test_plookup_n4() {
        let prime = BigUint::from_u32(31u32).unwrap();
        let f_values = vec![1, 1, 4, 8];
        let t_values = vec![1, 4, 8];

        let result = plookup_protocol(f_values, t_values, &prime).unwrap();

        assert!(result, "Plookup protocol failed for n = 4");
    }

    #[test]
    fn test_plookup_n9() {
        let prime = BigUint::from_u32(31u32).unwrap();
        let f_values = vec![1, 1, 1, 4, 8, 1, 4, 8, 1];
        let t_values = vec![1, 4, 8];

        let result = plookup_protocol(f_values, t_values, &prime).unwrap();

        assert!(result, "Plookup protocol failed for n = 9");
    }

    #[test]
    fn test_plookup_3() {
        let prime = BigUint::from_u32(31u32).unwrap();
        let f_values = vec![1, 2, 4, 8, 8];
        let t_values = vec![1, 2, 4, 8];

        let result = plookup_protocol(f_values, t_values, &prime).unwrap();

        assert!(result, "Plookup protocol failed for n = 4");
    }
}
