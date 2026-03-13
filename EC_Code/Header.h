#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdexcept>
#include <chrono>

#pragma warning(push)
#pragma warning(disable: 4146)  // unary minus on unsigned (GMP internal)
#include <gmpxx.h>
#pragma warning(pop)

// ======================== Modular Arithmetic ========================

// --- mod: повертає завжди невід'ємний результат ---

inline mpz_class mod(const mpz_class& value, const mpz_class& modulus) {
    mpz_class result = value % modulus;
    if (result < 0)
        result += modulus;
    return result;
}

inline long long mod(long long value, long long modulus) {
    long long result = value % modulus;
    if (result < 0)
        result += modulus;
    return result;
}

// --- modInverse: обернений елемент за модулем ---

// Спеціалізація для mpz_class (вбудована функція GMP)
inline mpz_class modInverse(const mpz_class& a, const mpz_class& p) {
    mpz_class result;
    if (mpz_invert(result.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t()) == 0)
        throw std::runtime_error("modInverse: inverse does not exist");
    return result;
}

// Спеціалізація для long long (розширений алгоритм Евкліда)
inline long long modInverse(long long a, long long p) {
    a = ((a % p) + p) % p;
    if (a == 0)
        throw std::runtime_error("modInverse: inverse of 0 does not exist");

    long long old_r = a, r = p;
    long long old_s = 1, s = 0;

    while (r != 0) {
        long long q = old_r / r;
        long long temp_r = r;
        r = old_r - q * r;
        old_r = temp_r;

        long long temp_s = s;
        s = old_s - q * s;
        old_s = temp_s;
    }

    if (old_r != 1)
        throw std::runtime_error("modInverse: inverse does not exist (gcd != 1)");

    return ((old_s % p) + p) % p;
}

// --- modPow: піднесення до степеня за модулем ---

// Спеціалізація для mpz_class (вбудована функція GMP)
inline mpz_class modPow(const mpz_class& base, const mpz_class& exp, const mpz_class& p) {
    mpz_class result;
    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    return result;
}

// Спеціалізація для long long (швидке піднесення до степеня)
inline long long modPow(long long base, long long exp, long long p) {
    base = ((base % p) + p) % p;
    long long result = 1;
    while (exp > 0) {
        if (exp % 2 == 1)
            result = (result * base) % p;
        base = (base * base) % p;
        exp /= 2;
    }
    return result;
}

// --- modSqrt: квадратний корінь за модулем простого p ---
// Повертає y таке що y^2 ≡ a (mod p), або -1 якщо кореня немає.
// Реалізація через алгоритм Тонеллі-Шенкса (працює для будь-якого непарного простого p).
// verbose = true -> покрокове виведення обчислень.

inline long long modSqrt(long long a, long long p, bool verbose = false) {
    a = mod(a, p);
    if (a == 0) return 0;

    // Перевірка: a — квадратичний лишок?
    // Символ Лежандра: (a/p) = a^((p-1)/2) mod p.  Результат: 1 = лишок, p-1 = нелишок.
    // Складність: O(log p) множень — швидко навіть для великих p.
    long long legendre = modPow(a, (p - 1) / 2, p);
    if (legendre != 1) {
        if (verbose) std::cout << "  [modSqrt] Legendre(" << a << "/" << p << ") = " << legendre << " != 1 -> not a QR" << std::endl;
        return -1;
    }

    // Випадок а): p ≡ 3 (mod 4) => y = a^((p+1)/4) mod p
    if (p % 4 == 3) {
        long long y = modPow(a, (p + 1) / 4, p);
        if (verbose)
            std::cout << "  [modSqrt] p = 4k+3 (k=" << (p - 3) / 4 << "): y = " << a
                      << "^(" << (p + 1) / 4 << ") = " << y << " (mod " << p << ")" << std::endl;
        return y;
    }

    // Випадок б): p ≡ 5 (mod 8) => формула Аткіна
    //   v = (2a)^((p-5)/8), i = 2a*v^2, y = a*v*(i-1)
    if (p % 8 == 5) {
        long long k = (p - 5) / 8;
        long long twoA = mod(2 * a, p);
        long long v = modPow(twoA, k, p);               // v = (2a)^k
        long long v2 = mod(v * v, p);                    // v^2
        long long i = mod(twoA * v2, p);                 // i = 2a*v^2
        long long y = mod(mod(a * v, p) * mod(i - 1, p), p);  // y = a*v*(i-1)
        if (verbose)
            std::cout << "  [modSqrt] p = 8k+5 (k=" << k << "): v=(2a)^k=" << v
                      << ", i=2a*v^2=" << i << ", y = a*v*(i-1) = " << y << " (mod " << p << ")" << std::endl;
        return y;
    }

    // Загальний випадок (p ≡ 1 mod 8): повний алгоритм Тонеллі-Шенкса
    if (verbose) std::cout << "  [modSqrt] p = 8k+1, using Tonelli-Shanks" << std::endl;

    // Розкладаємо p-1 = 2^s * q (q — непарне)
    long long s = 0, q = p - 1;
    while (q % 2 == 0) { q /= 2; s++; }

    // Шукаємо квадратичний нелишок z
    long long z = 2;
    while (modPow(z, (p - 1) / 2, p) != p - 1) z++;

    if (verbose)
        std::cout << "  p-1 = 2^" << s << " * " << q << ", non-residue z = " << z << std::endl;

    long long M = s;
    long long c = modPow(z, q, p);
    long long t = modPow(a, q, p);
    long long R = modPow(a, (q + 1) / 2, p);

    while (true) {
        if (t == 1) {
            if (verbose)
                std::cout << "  result y = " << R << " (mod " << p << ")" << std::endl;
            return R;
        }
        // Знаходимо найменше i таке що t^(2^i) ≡ 1
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

// Підключення класів кривої та точки
#include "EllipticCurve.h"
#include "EllipticCurvePoint.h"
