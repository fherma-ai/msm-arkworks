//! Multi-scalar multiplication on BLS12-381, through arkworks' own
//! variable-base routine.
//!
//! Two limb widths meet here, which is the one thing BN254 is spared. The base
//! field is 381 bits and takes six words; the scalar field is 255 and takes
//! four. The specification gives every tensor the same stride — six words to
//! a value — so a scalar arrives with its top two words zero and is narrowed
//! here rather than being read as a different shape.
use ark_bls12_381::{Fq, Fr, G1Affine, G1Projective};
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::{BigInt, PrimeField, Zero};

use crate::fherma::{Inputs, Outputs, Point, Tensor};

/// The `i`-th row of a limb tensor as a base-field element.
///
/// `from_bigint` refuses a value at or above the modulus. The generator
/// reduces before it writes, so a refusal here means the case is not the one
/// this specification describes — better said than folded away.
fn fq(data: &[u64], i: usize, limbs: usize) -> Fq {
    let mut words = [0u64; 6];
    words[..limbs].copy_from_slice(&data[i * limbs..(i + 1) * limbs]);
    Fq::from_bigint(BigInt::new(words)).expect("a base-field element below the modulus")
}

/// The `i`-th row as a scalar. The row is `limbs` words wide and the scalar
/// field is four; the words above it must be zero, and a case that sets one
/// is not this specification's.
fn fr(data: &[u64], i: usize, limbs: usize) -> Fr {
    let row = &data[i * limbs..(i + 1) * limbs];
    assert!(
        row[4.min(limbs)..].iter().all(|word| *word == 0),
        "a scalar wider than the scalar field",
    );
    let mut words = [0u64; 4];
    let take = limbs.min(4);
    words[..take].copy_from_slice(&row[..take]);
    Fr::from_bigint(BigInt::new(words)).expect("a scalar below the group order")
}

fn row(value: Fq, limbs: usize) -> Tensor<u64> {
    Tensor { shape: vec![limbs as i64], data: value.into_bigint().0[..limbs].to_vec() }
}

#[derive(Default)]
pub struct State {
    limbs: usize,
}

pub fn init(p: &Point) -> State {
    State { limbs: p.L as usize }
}

pub fn run(state: &mut State, inp: &Inputs) -> Outputs {
    let limbs = state.limbs;
    let n = inp.scalars.data.len() / limbs;

    let scalars: Vec<Fr> = (0..n).map(|i| fr(&inp.scalars.data, i, limbs)).collect();
    // `new_unchecked` because the generator takes every point from a square
    // root of the curve equation and then clears the cofactor, so each one is
    // in G1 by construction. Validating N points would be a second kernel in
    // front of this one.
    let bases: Vec<G1Affine> = (0..n)
        .map(|i| G1Affine::new_unchecked(fq(&inp.xs.data, i, limbs), fq(&inp.ys.data, i, limbs)))
        .collect();

    let sum = G1Projective::msm(&bases, &scalars)
        .expect("one scalar for every point")
        .into_affine();

    // The sum reaches infinity only by a collision nobody will see; the
    // specification writes it as (0, 0), which is not a point of this curve
    // and so cannot be mistaken for one.
    let (x, y) = sum.xy().unwrap_or((Fq::zero(), Fq::zero()));
    Outputs { rx: row(x, limbs), ry: row(y, limbs) }
}

pub fn free(state: State) {
    drop(state);
}
