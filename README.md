# Elliptic Curve Cryptography Library

C++ template library implementing elliptic curves in short Weierstrass form over prime fields,
with a focus on the Baby-JubJub curve.

---

## File Structure

```
Header.h              Modular arithmetic, polynomial arithmetic over F_p,
                      helper functions (Hasse bound, primality test, modSqrt)

EllipticCurve.h       Template class EllipticCurve<T> with curve parameters,
                      non-singularity check, bruteforce point enumeration,
                      and Schoof's algorithm for computing curve order

EllipticCurvePoint.h  Template class EllipticCurvePoint<T> with projective
                      and affine point arithmetic, scalar multiplication
                      (4 algorithms), point order computation

Main.cpp              Tests for simplified case (long long, small primes)
                      and Baby-JubJub (mpz_class, 254-bit prime),
                      including performance benchmarks
```

Both template classes support `T = long long` (small fields, fast, test) and `T = mpz_class`
(arbitrary precision via GMP, required for Baby-JubJub).

---

## Header.h -- Modular Arithmetic and Polynomials

### Modular arithmetic

| Function               | Description                                                                                                                                                                |
|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `mod(value, modulus)`  | Non-negative remainder. Always returns result in [0, modulus). Overloaded for `long long` and `mpz_class`.                                                                 |
| `modInverse(a, p)`     | Modular multiplicative inverse a^{-1} mod p. Uses GMP's `mpz_invert` for `mpz_class`, Extended Euclidean Algorithm for `long long`. Throws if inverse does not exist.      |
| `modPow(base, exp, p)` | Modular exponentiation base^exp mod p. Uses GMP's `mpz_powm` for `mpz_class`, binary square-and-multiply for `long long`.                                                  |
| `modSqrt(a, p)`        | Square root mod p (long long version). Returns y such that y^2 = a (mod p), or -1 if none exists. Uses Tonelli-Shanks with fast paths for p = 3 (mod 4) and p = 5 (mod 8). |
| `modSqrt(a, p)`        | Square root mod p (mpz_class version). Full Tonelli-Shanks for arbitrary-precision primes. Used for finding random curve points.                                           |

### Polynomial arithmetic over F_p

Polynomials are represented as `PolyT<T> = vector<T>`, where `poly[i]` = coefficient of x^i.
For example, 3x^2 + 5x + 1 is stored as {1, 5, 3}. The zero polynomial is {0} with degree -1
by convention (mathematically deg(0) = -infinity; we use -1 because it is sufficient to ensure
that the polynomial long division loop `while (deg(r) >= deg(b))` terminates when the
remainder becomes zero).

All polynomial operations work in the ring F_p[x] and are templated on T, so they support both
`long long` and `mpz_class` coefficients. This is critical for Schoof's algorithm, which needs
to compute x^p mod psi_l(x) where p can be a 254-bit prime.

| Function                        | Description                                                                                                                       |
|---------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `polyNorm(a)`                   | Remove leading zero coefficients. Ensures `poly.back() != 0` (except for zero polynomial).                                        |
| `polyDeg(a)`                    | Degree of polynomial. Returns -1 for zero polynomial (convention: deg(0) = -infinity).                                            |
| `polyAdd(a, b, p)`              | Coefficient-wise addition mod p.                                                                                                  |
| `polySub(a, b, p)`              | Coefficient-wise subtraction mod p.                                                                                               |
| `polyMul(a, b, p)`              | Naive schoolbook multiplication O(n*m). Each product is reduced mod p.                                                            |
| `polyScale(a, s, p)`            | Multiply polynomial by scalar s.                                                                                                  |
| `polyDivMod(a, b, p)`           | Long division: returns {quotient, remainder} such that a = q*b + r, deg(r) < deg(b).                                              |
| `polyMod(a, b, p)`              | Remainder of a mod b. Used for working in quotient ring F_p[x]/(b(x)).                                                            |
| `polyGcd(a, b, p)`              | Euclidean GCD algorithm. Result normalized to monic (leading coefficient = 1).                                                    |
| `polyPowMod(base, exp, mod, p)` | Binary exponentiation: base^exp mod modPoly in F_p[x]. Exponent has type T (not long long) to support 254-bit exponents like p^2. |
| `polySubst(poly, val, mod, p)`  | Composition: evaluate poly(val(x)) mod modPoly. Substitutes x -> val in poly with reduction after each step.                      |

### Helper functions

| Function          | Description                                                                                                 |
|-------------------|-------------------------------------------------------------------------------------------------------------|
| `hasseBound(p)`   | Computes 2*floor(sqrt(p)) + 1. Overloaded for `long long` (via `sqrt`) and `mpz_class` (via `mpz_sqrt`).    |
| `isSmallPrime(n)` | Deterministic primality test by trial division up to sqrt(n). Used for collecting small primes l in Schoof. |
| `fromLL<T>(v)`    | Convert `long long` to T without GMP's ambiguous `mpz_class(long long)` constructor.                        |
| `toLongLong(v)`   | Convert T to `long long` for values known to be small (e.g. p mod l).                                       |

---

## EllipticCurve.h -- Curve Class and Schoof's Algorithm

### Class: `EllipticCurve<T>`

Represents the curve y^2 = x^3 + ax + b over F_p with known (or to-be-computed) order n.

| Member    | Description                                                             |
|-----------|-------------------------------------------------------------------------|
| `a, b, p` | Curve coefficients and field prime                                      |
| `n`       | Curve order (number of points including O_E)                            |
| `verbose` | If true, print step-by-step computations (for long long demonstrations) |

### Methods

| Method                      | Description                                                                                                                                                    |
|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `isNonSingular()`           | Checks Delta = -16(4a^3 + 27b^2) != 0 mod p. Since -16 is invertible in F_p for p >= 3, only the core expression 4a^3 + 27b^2 is computed.                     |
| `print()`                   | Outputs the curve equation.                                                                                                                                    |
| `findAllPointsBruteforce()` | Iterates x from 0 to p-1, computes y^2 = f(x), attempts modSqrt. Stores order in n. Optionally prints a table of all points. Only for small p (T = long long). |
| `schoofOrder()`             | Computes curve order #E(F_p) = p + 1 - t using Schoof's algorithm. Returns the order as type T. Works for both long long and mpz_class.                        |

### Schoof's algorithm internals

The algorithm determines t mod l for small primes l, then combines results via CRT.

| Internal method                                | Description                                                                                                                                                                                                                                                  |
|------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `computeDivPolys(maxN)`                        | Computes division polynomials dp[0..maxN] using recurrence. dp[n] = psi_n for odd n, dp[n] = psi_n/(2y) for even n. Base cases dp[1]=1, dp[2]=1, dp[3], dp[4] are explicit formulas; dp[k] for k >= 5 uses the standard recurrence involving 16*f^2 factors. |
| `xCoordOfScalar(n, dp, fx, modP, num, den)`    | Computes x-coordinate of [n]P as a rational function num/den in F_p[x]/(modP). For odd n: num = x*dp[n]^2 - 4f*dp[n-1]*dp[n+1], den = dp[n]^2. For even n: the 4f factor moves to the denominator.                                                           |
| `yFactorOfScalar(n, dp, fx, modP, ynum, yden)` | Computes y-factor of [n]P: y([n]P) = y * ynum/yden. Common numerator for both parities: dp[n+2]*dp[n-1]^2 - dp[n-2]*dp[n+1]^2. Denominator: dp[n]^3 (odd n) or 16*f^2*dp[n]^3 (even n).                                                                      |

The main loop for each odd prime l:
1. Computes Frobenius x-coordinates: x^p and x^{p^2} mod psi_l
2. Computes Frobenius y-factors: f^{(p-1)/2} and f^{(p^2-1)/2} mod psi_l
3. Checks tau=0: whether pi^2(P) = -[q]P (or eigenvalue case pi^2(P) = +[q]P)
4. Iterates tau = 1..(l-1)/2: computes pi^2(P) + [q]P via chord formula in the polynomial ring, compares x-coordinate with [tau]*pi(P), then determines y-sign (tau vs l-tau)
5. CRT combines all t mod l values into t, then order = p + 1 - t

---

## EllipticCurvePoint.h -- Point Arithmetic

### Class: `EllipticCurvePoint<T>`

Represents a point on an elliptic curve in projective coordinates (X : Y : Z).
An affine point (x, y) is stored as (x, y, 1). The point at infinity O_E = (0 : 1 : 0).

### Constructors and conversion

| Method                             | Description                                                                    |
|------------------------------------|--------------------------------------------------------------------------------|
| `infinity(curve)`                  | Creates O_E = (0 : 1 : 0).                                                     |
| `fromAffine(x, y, curve)`          | Creates (x, y, 1) after validating the point lies on the curve. Throws if not. |
| `fromAffineUnchecked(x, y, curve)` | Creates (x, y, 1) without validation. For internal use by arithmetic results.  |
| `toAffine()`                       | Converts (X : Y : Z) to (X/Z, Y/Z) via modular inversion. Throws for O_E.      |
| `isInfinity()`                     | Returns true if Z == 0.                                                        |
| `isOnCurve()`                      | Checks Y^2*Z = X^3 + a*X*Z^2 + b*Z^3 mod p.                                    |
| `assertOnCurve(context)`           | Same check but throws with context message. Silent (no verbose output).        |

### Projective arithmetic

These formulas avoid modular inversion entirely, using only multiplications mod p.
The final result stays in projective form; conversion to affine (1 inversion) is done only when needed.

| Method          | Formula                                                                                                                        | Cost (approx.) |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------|----------------|
| `pointDouble()` | W = aZ^2 + 3X^2, S = YZ, B = XYS, H = W^2 - 8B, then X' = 2HS, Y' = W(4B-H) - 8Y^2S^2, Z' = 8S^3                               | 10M + 4S       |
| `pointAdd(Q)`   | U1=Y2*Z1, U2=Y1*Z2, V1=X2*Z1, V2=X1*Z2, U=U1-U2, V=V1-V2, W=Z1*Z2, A=U^2W-V^3-2V^2V2, then X3=VA, Y3=U(V^2V2-A)-V^3U2, Z3=V^3W | 12M + 2S       |

Handles special cases: P + O_E = P, P + (-P) = O_E, P + P delegates to pointDouble.

### Affine arithmetic

Each operation requires one modular inversion for the slope lambda.

| Method                | Formula                                                                            | Cost (approx.) |
|-----------------------|------------------------------------------------------------------------------------|----------------|
| `pointDoubleAffine()` | lambda = (3x^2 + a) / (2y), x' = lambda^2 - 2x, y' = lambda(x - x') - y            | 1I + 3M + 1S   |
| `pointAddAffine(Q)`   | lambda = (y2 - y1) / (x2 - x1), x' = lambda^2 - x1 - x2, y' = lambda(x1 - x') - y1 | 1I + 2M + 1S   |

If the point has Z != 1 (e.g. after projective operations), it is first converted to affine
via modInverse before applying the formula.

### Scalar multiplication

Four algorithms implementing kP, forming a 2x2 matrix of choices:

|            | Double-and-Add       | Montgomery Ladder              |
|------------|----------------------|--------------------------------|
| Projective | `scalarMul(k)`       | `scalarMulMontgomery(k)`       |
| Affine     | `scalarMulAffine(k)` | `scalarMulMontgomeryAffine(k)` |

Double-and-Add: iterates bits of k from LSB to MSB. On each step, doubles temp.
If bit is 1, adds temp to accumulator. Running time depends on Hamming weight of k (not constant-time).

Montgomery Ladder: iterates bits from MSB to LSB. Always performs exactly one addition
and one doubling per bit, regardless of the bit value. Constant-time execution, which is
important for cryptographic applications (protection against timing side-channel attacks).

### Other methods

| Method          | Description                                                                                                                                 |
|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| `pointOrder()`  | Finds the smallest k > 0 such that kP = O_E. By Lagrange's theorem, this divides the curve order n. Enumerates divisors of n from smallest. |
| `equals(other)` | Projective equality: checks X1*Z2 == X2*Z1 and Y1*Z2 == Y2*Z1 mod p.                                                                        |
| `print()`       | Outputs affine and projective coordinates.                                                                                                  |

---

## Main.cpp -- Tests and Benchmarks

### Simplified case (T = long long, p = 11)

Demonstrates all operations with verbose step-by-step output on the curve
y^2 = x^3 + 5x + 7 (mod 11):
- Bruteforce enumeration of all 16 points with sqrt computation example
- Schoof's algorithm verification (must match bruteforce)
- Point arithmetic: doubling, addition, inverse, identity element
- Scalar multiplication: 3P via all three algorithms, edge cases (0P, 1P)
- Affine verification: n*P = O_E using affine arithmetic
- Point order computation

### Baby-JubJub (T = mpz_class, p ~ 2^254)

- Derives Weierstrass parameters from Montgomery coefficient A = 168698 at runtime
- Converts the Edwards generator through Edwards -> Montgomery -> Weierstrass
- Verifies G + G = 2G, 5G via all algorithms
- Performance benchmark: 10 runs with random points, each generated using
  `gmp_randclass` seeded from OS entropy via `std::random_device`

### Performance: long long (p = 10^9 + 7)

Random points, scalar n*P where n ~ 10^9 (30 bits), average over 10 runs:

```
                       DoubleAndAdd    Montgomery
  Projective:             0.006 ms       0.007 ms
  Affine:                 0.009 ms       0.010 ms
```

### Performance: Baby-JubJub (p ~ 2^254)

Random points, scalar n*P where n ~ 2^254 (254 bits), average over 10 runs.

Results on Linux (GCC 13.3, GMP 6.3):

```
                         DoubleAndAdd    Montgomery
    Projective:              1.28 ms        2.56 ms
    Affine:                  1.00 ms        1.34 ms
```

Results on Windows (MSVC 2019, x64 Release, GMP via vcpkg):

```
                         DoubleAndAdd    Montgomery
    Projective:             21.25 ms       27.05 ms
    Affine:                  8.98 ms       14.40 ms
```

---

## Performance Analysis

### Why Montgomery is slower than Double-and-Add

Montgomery ladder always executes one addition + one doubling per bit, giving exactly
2 * 254 = 508 group operations for a 254-bit scalar. Double-and-Add performs 254 doublings +
~127 additions (on average, for random scalars) = ~381 operations. Montgomery pays ~33% more
operations for its constant-time guarantee. In a cryptographic setting this tradeoff is
worthwhile to prevent timing attacks, but for non-secret computations (like verifying n*P = O_E)
Double-and-Add is faster.

### Why affine is faster than projective for mpz_class

This result is counterintuitive, since the standard recommendation is to use projective
coordinates to avoid inversions. The explanation lies in GMP's implementation:

- GMP's `mpz_invert` uses a sub-quadratic half-GCD algorithm (Lehmer/Schoenhage).
  For 254-bit numbers, one inversion costs roughly the same as 5-8 multiplications.
- Projective doubling requires ~12 field multiplications per step but avoids inversion.
- Affine doubling requires ~4 multiplications + 1 inversion per step.
- Net cost per step: projective ~12M vs affine ~4M + 1I ~ 4M + 5-8M = 9-12M.

Since GMP's inversion is efficient at this operand size, affine comes out ahead.
Additionally, affine operations always produce Z = 1, and projective formulas
create more intermediate `mpz_class` temporary objects (W, S, B, H, etc.),
each requiring a heap allocation, adding constant-factor overhead.

### Why affine is slower for long long

For `long long`, our `modInverse` uses a hand-written Extended Euclidean Algorithm with a
loop of ~30 iterations (for 30-bit p). Although each iteration is a single CPU instruction,
the loop overhead and branch mispredictions make one inversion cost roughly 15-20
multiplications equivalent. This flips the balance:
projective ~12M vs affine ~4M + 1I ~ 4M + 15-20M = 19-24M, so projective wins.

### Note on X-only Montgomery ladder

An important optimization not implemented here: the Montgomery ladder can be performed using
only the x-coordinate, without tracking y at all. This is possible because the Montgomery
ladder maintains the invariant R1 - R0 = P throughout, and the addition formula for
x-coordinates on Montgomery curves requires only x(R0), x(R1), and x(P):

    x(R0 + R1) = x(P)^{-1} * (x(R0)*x(R1) - 1)^2 / (x(R0) - x(R1))^2

This "x-only" or "differential" addition avoids computing y entirely, saving roughly half
the multiplications per step. The y-coordinate can be recovered at the very end from x(P),
x(kP), and x((k+1)P) if needed. This approach is standard in real-world implementations
(e.g. X25519 for Curve25519) and would likely make the Montgomery ladder competitive with
or faster than Double-and-Add, especially for projective coordinates where the savings in
multiplications directly translate to fewer `mpz_class` temporaries.

---

## MSVC Warning Suppressions

| Warning | Location                      | Reason                                                                                                                                                                     |
|---------|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| C4146   | `#include <gmpxx.h>`          | GMP internally applies unary minus to unsigned types. Harmless.                                                                                                            |
| C26812  | `#include <gmpxx.h>`          | GMP uses unscoped `enum gmp_randalg_t`. Cannot be changed without modifying GMP headers.                                                                                   |
| C26451  | Polynomial arithmetic, Schoof | "int+int to size_t overflow". All polynomial indices are bounded by max(deg(psi_l)) = (l^2-1)/2. For Baby-JubJub, l < 100, so indices < 5000 -- far from INT_MAX (2*10^9). |
