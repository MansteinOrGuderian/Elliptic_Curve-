#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <algorithm>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4146)   // unary minus on unsigned (GMP internal)
#pragma warning(disable: 26812)  // unscoped enum in GMP headers (gmp_randalg_t)
#endif
#include <gmpxx.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

// ======================== Modular Arithmetic ========================

// --- mod: always returns a non-negative result in [0, modulus) ---

inline mpz_class mod(const mpz_class& value, const mpz_class& modulus) {
    mpz_class result = value % modulus;
    if (result < 0) result += modulus;
    return result;
}

inline long long mod(long long value, long long modulus) {
    long long result = value % modulus;
    if (result < 0) result += modulus;
    return result;
}

// --- modInverse: modular multiplicative inverse (a^{-1} mod p) ---
// Uses GMP's built-in mpz_invert for mpz_class,
// and the Extended Euclidean Algorithm for long long.

inline mpz_class modInverse(const mpz_class& a, const mpz_class& p) {
    mpz_class result;
    if (mpz_invert(result.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t()) == 0)
        throw std::runtime_error("modInverse: inverse does not exist");
    return result;
}

inline long long modInverse(long long a, long long p) {
    a = ((a % p) + p) % p;
    if (a == 0)
        throw std::runtime_error("modInverse: inverse of 0 does not exist");
    long long old_r = a, r = p;
    long long old_s = 1, s = 0;
    while (r != 0) {
        long long q = old_r / r;
        long long temp_r = r; r = old_r - q * r; old_r = temp_r;
        long long temp_s = s; s = old_s - q * s; old_s = temp_s;
    }
    if (old_r != 1)
        throw std::runtime_error("modInverse: inverse does not exist (gcd != 1)");
    return ((old_s % p) + p) % p;
}

// --- modPow: modular exponentiation (base^exp mod p) ---
// Uses GMP's built-in mpz_powm for mpz_class,
// and binary exponentiation (square-and-multiply) for long long.

inline mpz_class modPow(const mpz_class& base, const mpz_class& exp, const mpz_class& p) {
    mpz_class result;
    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    return result;
}

inline long long modPow(long long base, long long exp, long long p) {
    base = ((base % p) + p) % p;
    long long result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) result = (result * base) % p;
        base = (base * base) % p;
        exp /= 2;
    }
    return result;
}

// --- modSqrt: modular square root over a prime field ---
// Returns y such that y^2 = a (mod p), or -1 if no root exists.
// Uses Tonelli-Shanks algorithm with fast paths for p = 3 (mod 4) and p = 5 (mod 8).

inline long long modSqrt(long long a, long long p, bool verbose = false) {
    a = mod(a, p);
    if (a == 0) return 0;

    // Legendre symbol: a^{(p-1)/2} mod p.  Returns 1 if QR, p-1 if NQR.
    long long legendre = modPow(a, (p - 1) / 2, p);
    if (legendre != 1) {
        if (verbose) std::cout << "  [modSqrt] Legendre(" << a << "/" << p << ") = " << legendre << " != 1 -> not a QR" << std::endl;
        return -1;
    }

    // Case a) p = 3 (mod 4): y = a^{(p+1)/4}
    if (p % 4 == 3) {
        long long y = modPow(a, (p + 1) / 4, p);
        if (verbose)
            std::cout << "  [modSqrt] p = 4k+3: y = " << a << "^(" << (p + 1) / 4 << ") = " << y << " (mod " << p << ")" << std::endl;
        return y;
    }

    // Case b) p = 5 (mod 8): Atkin's formula
    if (p % 8 == 5) {
        long long k = (p - 5) / 8;
        long long twoA = mod(2 * a, p);
        long long v = modPow(twoA, k, p);
        long long v2 = mod(v * v, p);
        long long i = mod(twoA * v2, p);
        long long y = mod(mod(a * v, p) * mod(i - 1, p), p);
        if (verbose)
            std::cout << "  [modSqrt] p = 8k+5: v=(2a)^k=" << v << ", i=2a*v^2=" << i << ", y = a*v*(i-1) = " << y << std::endl;
        return y;
    }

    // General case: full Tonelli-Shanks (p = 1 mod 8)
    if (verbose) std::cout << "  [modSqrt] p = 8k+1, using Tonelli-Shanks" << std::endl;

    // Factor p-1 = 2^s * q (q odd)
    long long s = 0, q = p - 1;
    while (q % 2 == 0) { q /= 2; s++; }
    // Find a quadratic non-residue z
    long long z = 2;
    while (modPow(z, (p - 1) / 2, p) != p - 1) z++;
    if (verbose) std::cout << "  p-1 = 2^" << s << " * " << q << ", non-residue z = " << z << std::endl;

    long long M = s;
    long long c = modPow(z, q, p);
    long long t = modPow(a, q, p);
    long long R = modPow(a, (q + 1) / 2, p);

    while (true) {
        if (t == 1) {
            if (verbose) std::cout << "  result y = " << R << " (mod " << p << ")" << std::endl;
            return R;
        }
        // Find smallest i such that t^{2^i} = 1
        long long i = 0;
        long long temp = t;
        while (temp != 1) { temp = mod(temp * temp, p); i++; }
        long long b = c;
        for (long long j = 0; j < M - i - 1; j++) b = mod(b * b, p);
        M = i;
        c = mod(b * b, p);
        t = mod(t * c, p);
        R = mod(R * b, p);
    }
}

// ======================== Templated Polynomial Arithmetic over F_p ========================
//
// Polynomials are represented as vector<T>, where poly[i] = coefficient of x^i.
// Example: 3x^2 + 5x + 1 is stored as {1, 5, 3}.
// Zero polynomial: {0}.  Degree of zero polynomial = -1 (by convention).
//
// All operations are performed modulo a prime p, i.e. we work in the ring F_p[x].
//
// These functions are used by Schoof's algorithm, which operates in the quotient ring
// F_p[x]/(psi_l(x)), where psi_l is the l-th division polynomial.

// MSVC C26451: "Arithmetic overflow: int+int -> size_t". This is a false positive here.
// All polynomial indices are bounded by max(deg(psi_l)) = (l^2-1)/2.
// For Schoof on Baby-JubJub (p ~ 2^254), l < 100, so indices < 5000 — far from INT_MAX.
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 26451)
#endif

template<typename T>
using PolyT = std::vector<T>;

// Backward compatibility alias: Poly = PolyT<long long>
using Poly = PolyT<long long>;

// --- polyNorm: normalize polynomial ---
// Removes leading zero coefficients: {1, 0, 3, 0, 0} -> {1, 0, 3}.
// Guarantees poly.back() != 0 (except for zero polynomial {0}).
// Without normalization, degree computation would be incorrect.
template<typename T>
PolyT<T> polyNorm(PolyT<T> a) {
    while (a.size() > 1 && a.back() == T(0)) a.pop_back();
    return a;
}

// --- polyDeg: polynomial degree ---
// Returns the highest power of x with a nonzero coefficient.
// Zero polynomial {0} returns -1 (convention: deg(0) = -inf).
// Critical for polyDivMod: the division loop terminates when deg(r) < deg(b).
template<typename T>
int polyDeg(const PolyT<T>& a) {
    if (a.size() == 1 && a[0] == T(0)) return -1;
    return static_cast<int>(a.size()) - 1;
}

// --- polyAdd: polynomial addition ---
// (a0 + a1*x + ...) + (b0 + b1*x + ...) = (a0+b0) + (a1+b1)*x + ...
// Coefficient-wise addition modulo p. Result is normalized.
template<typename T>
PolyT<T> polyAdd(const PolyT<T>& a, const PolyT<T>& b, const T& p) {
    PolyT<T> res(std::max(a.size(), b.size()), T(0));
    for (size_t i = 0; i < a.size(); i++) res[i] = mod(res[i] + a[i], p);
    for (size_t i = 0; i < b.size(); i++) res[i] = mod(res[i] + b[i], p);
    return polyNorm(res);
}

// --- polySub: polynomial subtraction ---
// a(x) - b(x), coefficient-wise modulo p.
template<typename T>
PolyT<T> polySub(const PolyT<T>& a, const PolyT<T>& b, const T& p) {
    PolyT<T> res(std::max(a.size(), b.size()), T(0));
    for (size_t i = 0; i < a.size(); i++) res[i] = mod(res[i] + a[i], p);
    for (size_t i = 0; i < b.size(); i++) res[i] = mod(res[i] - b[i], p);
    return polyNorm(res);
}

// --- polyMul: polynomial multiplication ---
// Naive schoolbook multiplication O(n*m):
//   (sum a_i x^i) * (sum b_j x^j) = sum (a_i * b_j) * x^{i+j}
// Degree of result = deg(a) + deg(b).
// Each term is reduced mod p to prevent overflow.
template<typename T>
PolyT<T> polyMul(const PolyT<T>& a, const PolyT<T>& b, const T& p) {
    if (a.empty() || b.empty()) return {T(0)};
    PolyT<T> res(a.size() + b.size() - 1, T(0));
    for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < b.size(); j++)
            res[i + j] = mod(res[i + j] + mod(a[i] * b[j], p), p);
    return polyNorm(res);
}

// --- polyScale: multiply polynomial by a scalar ---
// s * (a0 + a1*x + a2*x^2 + ...) = (s*a0) + (s*a1)*x + (s*a2)*x^2 + ...
template<typename T>
PolyT<T> polyScale(const PolyT<T>& a, const T& s, const T& p) {
    PolyT<T> res(a.size());
    for (size_t i = 0; i < a.size(); i++) res[i] = mod(a[i] * s, p);
    return polyNorm(res);
}

// --- polyDivMod: polynomial long division with remainder ---
//
// For a(x) and b(x) != 0, finds q(x) and r(x) such that:
//     a(x) = q(x) * b(x) + r(x),   deg(r) < deg(b)
//
// Algorithm (analogous to integer long division):
//   1. While deg(r) >= deg(b):
//      - c = leading_coeff(r) / leading_coeff(b)
//      - r -= c * x^shift * b(x),  where shift = deg(r) - deg(b)
//   2. Returns: {quotient, remainder}.
template<typename T>
std::pair<PolyT<T>, PolyT<T>> polyDivMod(const PolyT<T>& a, const PolyT<T>& b, const T& p) {
    if (polyDeg(b) < 0 || (b.size() == 1 && b[0] == T(0)))
        throw std::runtime_error("polyDivMod: division by zero polynomial");

    PolyT<T> r = a;
    int degB = polyDeg(b);
    T leadInv = modInverse(b[degB], p);
    PolyT<T> q(std::max((int)a.size() - degB, 1), T(0));

    while (polyDeg(polyNorm(r)) >= degB) {
        r = polyNorm(r);
        int degR = polyDeg(r);
        T coeff = mod(r[degR] * leadInv, p);
        int shift = degR - degB;
        q[shift] = coeff;
        for (int i = 0; i <= degB; i++)
            r[i + shift] = mod(r[i + shift] - mod(coeff * b[i], p), p);
    }
    return {polyNorm(q), polyNorm(r)};
}

// --- polyMod: polynomial remainder ---
// a(x) mod b(x) — used for working in the quotient ring F_p[x]/(b(x)).
template<typename T>
PolyT<T> polyMod(const PolyT<T>& a, const PolyT<T>& b, const T& p) {
    return polyDivMod(a, b, p).second;
}

// --- polyGcd: greatest common divisor (Euclidean algorithm) ---
// gcd(a, b) = gcd(b, a mod b), gcd(a, 0) = a.
// Result is normalized to monic (leading coefficient = 1).
//
// Used in Schoof: for l=2, we check gcd(x^p - x, f(x)).
// If degree > 0, f has a root in F_p, hence t = 0 (mod 2).
template<typename T>
PolyT<T> polyGcd(PolyT<T> a, PolyT<T> b, const T& p) {
    while (!(b.size() == 1 && b[0] == T(0))) {
        PolyT<T> r = polyMod(a, b, p);
        a = b;
        b = r;
    }
    if (!a.empty() && a.back() != T(0)) {
        T inv = modInverse(a.back(), p);
        for (auto& c : a) c = mod(c * inv, p);
    }
    return polyNorm(a);
}

// --- polyPowMod: binary exponentiation for polynomials ---
// Computes base(x)^exp mod modPoly(x) in F_p[x].
//
// The exponent has type T (not long long!) because Schoof needs to compute
// x^p mod psi_l(x) and f(x)^{(p^2-1)/2} mod psi_l(x), where p can be 254-bit.
// Complexity: O(log(exp)) polynomial multiplications.
template<typename T>
PolyT<T> polyPowMod(PolyT<T> base, T exp, const PolyT<T>& modPoly, const T& p) {
    PolyT<T> result = {T(1)};
    base = polyMod(base, modPoly, p);
    while (exp > T(0)) {
        if (exp % T(2) == T(1))
            result = polyMod(polyMul(result, base, p), modPoly, p);
        base = polyMod(polyMul(base, base, p), modPoly, p);
        exp /= T(2);
    }
    return result;
}

// --- polySubst: polynomial composition ---
// Computes poly(val(x)) mod modPoly(x) in F_p[x].
// Substitutes x -> val(x) in poly, reducing mod modPoly after each step.
//
// Example: poly = 3x^2 + 2x + 1, val = x^3 -> 3x^6 + 2x^3 + 1.
//
// Used in Schoof: substitute x^p (Frobenius) into division polynomials dp[tau]
// to compute dp[tau](x^p) — the x-coordinate of [tau]*pi(P) in F_p[x]/(psi_l).
template<typename T>
PolyT<T> polySubst(const PolyT<T>& poly, const PolyT<T>& val,
                    const PolyT<T>& modPoly, const T& p) {
    PolyT<T> result = {T(0)};
    PolyT<T> val_pow = {T(1)};
    for (size_t i = 0; i < poly.size(); i++) {
        if (poly[i] != T(0))
            result = polyAdd(result, polyScale(val_pow, poly[i], p), p);
        if (!modPoly.empty()) result = polyMod(result, modPoly, p);
        if (i + 1 < poly.size()) {
            val_pow = polyMul(val_pow, val, p);
            if (!modPoly.empty()) val_pow = polyMod(val_pow, modPoly, p);
        }
    }
    return result;
}

// ======================== Helper Functions ========================

// Hasse bound: |t| <= 2*floor(sqrt(p)) + 1, where #E = p + 1 - t.
inline long long hasseBound(long long p) {
    return (long long)(2.0 * std::sqrt((double)p)) + 1;
}
inline mpz_class hasseBound(const mpz_class& p) {
    mpz_class s;
    mpz_sqrt(s.get_mpz_t(), p.get_mpz_t());
    return 2 * s + 1;
}

// Deterministic primality test for small numbers (trial division up to sqrt(n)).
inline bool isSmallPrime(long long n) {
    if (n < 2) return false;
    if (n < 4) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (long long d = 5; d * d <= n; d += 6)
        if (n % d == 0 || n % (d + 2) == 0) return false;
    return true;
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Include curve and point class definitions
#include "EllipticCurve.h"
#include "EllipticCurvePoint.h"
