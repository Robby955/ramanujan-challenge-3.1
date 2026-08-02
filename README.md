# The $7_2$ regulator integral is $4\pi^2/85$

## A reproducible computer-assisted proof of the $7_2$ regulator identity

Status: **submitted; organizer evaluation pending**.

The Ramanujan Challenge closed on August 1, 2026. One of the open problems it
posed, Problem 3.1, asked for the exact value of a regulator integral associated
with the $7_2$ knot. This repository contains a computer-assisted proof that the
value is $4\pi^2/85$.

Rob Sneiderman submitted the proof on 2026-07-12 at 02:38:40 UTC. The
submission is recorded by the organizers under SHA-256 receipt
`87985f32...9f06660`.

## Result

For the distinguished real branch of the `7_2` knot A-polynomial specified in
[Problem 3.1](https://www.ramanujanmachine.com/ramanujan-challenge/), the
manuscript proves

```text
Integral_alpha^beta (log(x) dy/y - log(y) dx/x) = 4*pi^2/85.
```

The challenge introduced this formula as an open conjecture. Khoi's 2008 paper
develops the general Godbillon-Vey/A-polynomial strategy and evaluates related
paths for the figure-eight and `5_2` knots, but not this `7_2` identity. The
claimed contribution here is the explicit `7_2` proof and reproducible
certificate, not a new general regulator method.

## Proof spine

1. An exact rational formula lifts the full A-polynomial arc to a continuous,
   irreducible Riley representation path.
2. The endpoint characters factor through the `(-1,2)` and `(-1,1)` Dehn
   fillings.
3. The four-tetrahedron triangulation, gluing rows, positive algebraic shapes,
   and flattenings are reconstructed without SnapPy's private
   `dev.extended_ptolemy` interface.
4. Reciprocal Galois symmetry descends twice each endpoint Bloch class to a
   totally real field. Borel finiteness and the cited regulator normalization
   place the normalized difference on a lattice of spacing `1/17821440`.
5. A 1000-bit Arb enclosure of diameter about `5e-300` contains `4/85` and
   rigorously excludes both neighboring lattice points.
6. The independently reconstructed Dehn wedge gives
   `dS = log(M) dlog(L) - log(L) dlog(M)` with the required sign. Positivity
   and `0<I<pi^2/2` select the real lift.

## Main artifacts

- `paper/manuscript.pdf`: public-facing manuscript.
- `paper/manuscript.tex`: manuscript source.
- `verification/LEDGER.md`: claim-by-claim dependency ledger.
- `verification/run_all.sh`: complete replay gate.
- `verification/exactness_closure.py`: final lattice and certified
  interval selection.
- `verification/byhand_flattening.py`: triangulation, gluing,
  shape, flattening, and regulator reconstruction.
- `submission/original-submission.pdf`: original seven-page challenge
  submission, retained as the submission record.

## Reproduce

Requirements:

- Python 3 with `sympy` and `mpmath`;
- SageMath 10.9;
- SnapPy 3.3.2 installed in the same Python interpreter used by Sage.

Install the plain-Python dependencies from the repository root:

```bash
python3 -m pip install -r requirements.txt
```

If Sage's Python does not already contain SnapPy:

```bash
export SAGE=/path/to/sage
"$SAGE" -python -m pip install snappy==3.3.2
```

Then run from the repository root:

```bash
verification/run_all.sh
```

The final script prints
`CLOSED-FULLY (modulo cited published theorems)` only after the exact field,
torsion-bound, Arb-enclosure, and adjacent-lattice-point assertions pass.
Borel's theorem and the Zickert-Neumann normalization are cited inputs rather
than re-proved results.

## Current audit status

The complete independent suite was replayed on 2026-08-01. The exact Riley and
Ptolemy identities, endpoint uniqueness, Dehn-wedge sign, 60-digit direct
quadrature, by-hand regulator reconstruction, and 1000-bit lattice selection
all passed. During the replay, a recursive Sage launcher failure was found and
repaired; the scripts now fail with an installation instruction instead of
re-executing indefinitely when SnapPy is missing from Sage's Python.

The audit establishes that the proof is reproducible from the released
artifacts. Organizer evaluation and further independent mathematical review
remain pending.

## License

The manuscript and written material are licensed under CC BY 4.0. The
verification code is licensed under the MIT License. See `LICENSE`,
`LICENSE-PAPER`, and `LICENSE-CODE`.

## Tooling disclosure

The project was developed by Rob Sneiderman with extensive assistance from
OpenAI Codex and other frontier language models for symbolic exploration,
implementation, adversarial review, and exposition. Model output is not treated
as proof; the mathematical claim rests on the manuscript, cited theorems, and
reproducible exact and certified computations.
