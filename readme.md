# Lookup Table

This is a tiny Rust implementation written to understand primary lookup tables in ZKP. **This implementation is not secure and is intended solely for my practice purposes.** It doesn't support the commitment scheme, and focus on the operations of polynomials.

## Features

- **Plookup**: Implements the plookup protocol using univariate polynomials.
- **Logup**: Implements the logup protocol based on the sumcheck protocol. It doesn't support for the GKR protocol yet.

## Components

- `field.rs`: Implements basic arithmetic operations on finite fields.
- `univariate_polynomial.rs`: Provides functionality for creating and manipulating univariate polynomials.
- `plookup.rs`: Implements the plookup protocol using univariate polynomials.
- `multivariate_polynomial.rs`: Handles multivariate polynomials required for the sumcheck protocol.
- `sumcheck.rs`: Implements the sumcheck protocol, a fundamental component of the logup protocol.
- `logup.rs`: Implements the logup protocol based on the sumcheck protocol. 

## Reference

[plookup: A simplified polynomial protocol for lookup tables](https://eprint.iacr.org/2020/315.pdf)

[Multivariate lookups based on logarithmic derivatives](https://eprint.iacr.org/2022/1530.pdf)

[Improving logarithmic derivative lookups using
GKR](https://eprint.iacr.org/2023/1284.pdf)

[pylookup](https://github.com/NOOMA-42/pylookup)