//! Multi-scalar multiplication, through arkworks' own variable-base routine.
use ark_bn254::{Fq, Fr, G1Affine, G1Projective};
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::{BigInt, PrimeField, Zero};

use crate::fherma::{Inputs, Outputs, Point, Tensor};

/// The `i`-th row of a limb tensor as an element of the field `F`.
///
/// `from_bigint` refuses a value at or above the modulus. The generator
/// reduces before it writes, so a refusal here means the case is not the one
/// this specification describes — better said than folded away.
fn fe<F: PrimeField<BigInt = BigInt<4>>>(data: &[u64], i: usize, limbs: usize) -> F {
    let mut words = [0u64; 4];
    words[..limbs].copy_from_slice(&data[i * limbs..(i + 1) * limbs]);
    F::from_bigint(BigInt::new(words)).expect("a field element below the modulus")
}

fn row<F: PrimeField<BigInt = BigInt<4>>>(value: F, limbs: usize) -> Tensor<u64> {
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

    let scalars: Vec<Fr> = (0..n).map(|i| fe(&inp.scalars.data, i, limbs)).collect();
    // `new_unchecked` because the generator takes every point from a square
    // root of the curve equation, so they are on the curve by construction and
    // BN254 has cofactor one — there is no subgroup left to check. Validating
    // N points would be a second kernel in front of this one.
    let bases: Vec<G1Affine> = (0..n)
        .map(|i| {
            G1Affine::new_unchecked(
                fe::<Fq>(&inp.xs.data, i, limbs),
                fe::<Fq>(&inp.ys.data, i, limbs),
            )
        })
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
