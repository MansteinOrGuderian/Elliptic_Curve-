#pragma once

// Elliptic curve point in projective coordinates (X, Y, Z)
// Affine point (x, y) corresponds to projective (x, y, 1)
// Point at infinity O_E = (0, 1, 0)

template <typename T>
class EllipticCurvePoint {
public:
    T X, Y, Z;
    const EllipticCurve<T>* curve;  // pointer to the curve this point belongs to

    // --- Constructors ---

    // Projective point constructor
    EllipticCurvePoint(const T& X, const T& Y, const T& Z, const EllipticCurve<T>* curve)
        : X(X), Y(Y), Z(Z), curve(curve) {}

    // Create the point at infinity
    static EllipticCurvePoint infinity(const EllipticCurve<T>* curve) {
        return EllipticCurvePoint(T(0), T(1), T(0), curve);
    }

    // Create a point from affine coordinates (x, y) -> (x, y, 1)
    // Validates that the point lies on the curve; throws if not.
    static EllipticCurvePoint fromAffine(const T& x, const T& y, const EllipticCurve<T>* curve) {
        EllipticCurvePoint pt(x, y, T(1), curve);
        if (!pt.isOnCurve())
            throw std::runtime_error("fromAffine: point (" + toString(x) + ", " + toString(y) + ") is NOT on the curve");
        return pt;
    }

    // Create a point WITHOUT validation (for internal use —
    // arithmetic results are guaranteed on the curve if the input was correct)
    static EllipticCurvePoint fromAffineUnchecked(const T& x, const T& y, const EllipticCurve<T>* curve) {
        return EllipticCurvePoint(x, y, T(1), curve);
    }

    // --- Checks ---

    bool isInfinity() const {
        return Z == 0;
    }

    // Throwing check — used before arithmetic operations.
    // Silent: does not print verbose output (internal check).
    void assertOnCurve(const std::string& context) const {
        if (isInfinity()) return;  // O_E is always on the curve

        const T& p = curve->p;
        const T& a = curve->a;
        const T& b = curve->b;

        T Y2 = mod(Y * Y, p);
        T lhs = mod(Y2 * Z, p);
        T X2 = mod(X * X, p);
        T X3 = mod(X2 * X, p);
        T Z2 = mod(Z * Z, p);
        T Z3 = mod(Z2 * Z, p);
        T rhs = mod(X3 + mod(mod(a * X, p) * Z2, p) + mod(b * Z3, p), p);

        if (lhs != rhs)
            throw std::runtime_error(context + ": point is NOT on the curve");
    }

    // Check whether the point lies on the curve: Y^2*Z = X^3 + a*X*Z^2 + b*Z^3 (mod p)
    bool isOnCurve() const {
        if (isInfinity()) return true;

        const T& p = curve->p;
        const T& a = curve->a;
        const T& b = curve->b;

        // mod() after every multiplication to prevent long long overflow
        T Y2 = mod(Y * Y, p);
        T lhs = mod(Y2 * Z, p);                                // Y^2 * Z

        T X2 = mod(X * X, p);
        T X3 = mod(X2 * X, p);                                 // X^3
        T Z2 = mod(Z * Z, p);
        T Z3 = mod(Z2 * Z, p);                                 // Z^3
        T aXZ2 = mod(mod(a * X, p) * Z2, p);                   // a*X*Z^2
        T bZ3  = mod(b * Z3, p);                                // b*Z^3
        T rhs  = mod(X3 + aXZ2 + bZ3, p);                      // X^3 + a*X*Z^2 + b*Z^3

        if (curve->verbose) {
            std::cout << "[isOnCurve] (" << toString(X) << ", " << toString(Y) << ", " << toString(Z) << ")" << std::endl;
            std::cout << "  Left:  Y^2*Z = " << toString(Y) << "^2 * " << toString(Z) << " = " << toString(lhs) << " (mod " << toString(p) << ")" << std::endl;
            std::cout << "  Right: X^3 + a*X*Z^2 + b*Z^3 = " << toString(rhs) << " (mod " << toString(p) << ")" << std::endl;
            std::cout << "  -> " << (lhs == rhs ? "ON curve" : "NOT on curve") << std::endl;
        }

        return lhs == rhs;
    }

    // --- Conversion ---

    // Projective -> affine: (X, Y, Z) -> (X/Z, Y/Z)
    // Returns a pair (x, y). Throws for O_E (no affine representation).
    std::pair<T, T> toAffine() const {
        if (isInfinity())
            throw std::runtime_error("toAffine: point at infinity has no affine coordinates");

        const T& p = curve->p;
        T zInv = modInverse(Z, p);
        T x = mod(X * zInv, p);
        T y = mod(Y * zInv, p);

        if (curve->verbose) {
            std::cout << "[toAffine] (" << toString(X) << " : " << toString(Y) << " : " << toString(Z) << ") -> ("
                      << toString(x) << ", " << toString(y) << ")" << std::endl;
        }

        return { x, y };
    }

    // --- Arithmetic ---

    // Point doubling in projective coordinates.
    // Formulas: W = a*Z^2 + 3*X^2, S = Y*Z, B = X*Y*S,
    //   H = W^2 - 8*B, X' = 2*H*S, Y' = W*(4*B - H) - 8*Y^2*S^2, Z' = 8*S^3
    EllipticCurvePoint pointDouble() const {
        if (isInfinity()) {
            if (curve->verbose)
                std::cout << "[pointDouble] O_E doubled -> O_E" << std::endl;
            return infinity(curve);
        }

        assertOnCurve("pointDouble");

        const T& p = curve->p;
        const T& a = curve->a;

        // Point of order 2: Y == 0 => 2P = O_E
        if (mod(Y, p) == 0) {
            if (curve->verbose)
                std::cout << "[pointDouble] P = (" << toString(X) << ", " << toString(Y) << ", " << toString(Z)
                          << "), Y == 0 -> point of order 2 -> O_E" << std::endl;
            return infinity(curve);
        }

        // Projective doubling formulas.
        // mod() after every multiplication to prevent long long overflow.
        T W  = mod(mod(a * mod(Z * Z, p), p) + mod(3 * mod(X * X, p), p), p);
        T S  = mod(Y * Z, p);
        T B  = mod(mod(X * Y, p) * S, p);
        T H  = mod(mod(W * W, p) - mod(8 * B, p), p);
        T Xr = mod(mod(2 * H, p) * S, p);
        T Yr = mod(mod(W * mod(4 * B - H, p), p)
                  - mod(mod(8 * mod(Y * Y, p), p) * mod(S * S, p), p), p);
        T Zr = mod(8 * mod(mod(S * S, p) * S, p), p);

        if (curve->verbose) {
            std::string mp = " (mod " + toString(p) + ")";
            std::cout << "[pointDouble] P = (" << toString(X) << ", " << toString(Y) << ", " << toString(Z) << ")" << std::endl;
            std::cout << "  W = a*Z^2 + 3*X^2 = " << toString(a) << "*" << toString(Z) << "^2 + 3*" << toString(X) << "^2 = " << toString(W) << mp << std::endl;
            std::cout << "  S = Y*Z = " << toString(Y) << "*" << toString(Z) << " = " << toString(S) << mp << std::endl;
            std::cout << "  B = X*Y*S = " << toString(X) << "*" << toString(Y) << "*" << toString(S) << " = " << toString(B) << mp << std::endl;
            std::cout << "  H = W^2 - 8*B = " << toString(W) << "^2 - 8*" << toString(B) << " = " << toString(H) << mp << std::endl;
            std::cout << "  X' = 2*H*S = 2*" << toString(H) << "*" << toString(S) << " = " << toString(Xr) << mp << std::endl;
            std::cout << "  Y' = W*(4B-H) - 8*Y^2*S^2 = " << toString(Yr) << mp << std::endl;
            std::cout << "  Z' = 8*S^3 = 8*" << toString(S) << "^3 = " << toString(Zr) << mp << std::endl;
        }

        return EllipticCurvePoint(Xr, Yr, Zr, curve);
    }

    // Point addition in projective coordinates.
    // Formulas: U1, U2, V1, V2 -> U, V, W, A -> X3, Y3, Z3
    EllipticCurvePoint pointAdd(const EllipticCurvePoint& other) const {
        const T& p = curve->p;

        // P + O_E = P,  O_E + Q = Q
        if (isInfinity()) {
            if (curve->verbose)
                std::cout << "[pointAdd] P = O_E -> result = Q" << std::endl;
            return other;
        }
        if (other.isInfinity()) {
            if (curve->verbose)
                std::cout << "[pointAdd] Q = O_E -> result = P" << std::endl;
            return *this;
        }

        assertOnCurve("pointAdd (P)");
        other.assertOnCurve("pointAdd (Q)");

        // U1 = Y2*Z1, U2 = Y1*Z2 (compare Y-coordinates in common scale)
        T U1 = mod(other.Y * Z, p);
        T U2 = mod(Y * other.Z, p);
        // V1 = X2*Z1, V2 = X1*Z2 (compare X-coordinates in common scale)
        T V1 = mod(other.X * Z, p);
        T V2 = mod(X * other.Z, p);

        if (V1 == V2) {
            // Same X-coordinate
            if (U1 != U2) {
                // P and Q are inverses: P + (-P) = O_E
                if (curve->verbose)
                    std::cout << "[pointAdd] P = (" << toString(X) << ", " << toString(Y) << ", " << toString(Z)
                              << "), Q = (" << toString(other.X) << ", " << toString(other.Y) << ", " << toString(other.Z)
                              << ") -> inverse points -> O_E" << std::endl;
                return infinity(curve);
            } else {
                // P == Q: delegate to point doubling
                if (curve->verbose)
                    std::cout << "[pointAdd] P == Q -> calling pointDouble" << std::endl;
                return pointDouble();
            }
        }

        // General case: P != Q, P != -Q
        T U  = mod(U1 - U2, p);
        T V  = mod(V1 - V2, p);
        T W  = mod(Z * other.Z, p);
        T Vsq = mod(V * V, p);
        T Vcb = mod(Vsq * V, p);
        T Usq = mod(U * U, p);
        T VsqV2 = mod(Vsq * V2, p);
        T A  = mod(mod(Usq * W, p) - Vcb - mod(2 * VsqV2, p), p);
        T X3 = mod(V * A, p);
        T Y3 = mod(mod(U * mod(VsqV2 - A, p), p) - mod(Vcb * U2, p), p);
        T Z3 = mod(Vcb * W, p);

        if (curve->verbose) {
            std::string mp = " (mod " + toString(p) + ")";
            std::cout << "[pointAdd] P = (" << toString(X) << ", " << toString(Y) << ", " << toString(Z)
                      << "), Q = (" << toString(other.X) << ", " << toString(other.Y) << ", " << toString(other.Z) << ")" << std::endl;
            std::cout << "  U1 = Y2*Z1 = " << toString(other.Y) << "*" << toString(Z) << " = " << toString(U1) << mp << std::endl;
            std::cout << "  U2 = Y1*Z2 = " << toString(Y) << "*" << toString(other.Z) << " = " << toString(U2) << mp << std::endl;
            std::cout << "  V1 = X2*Z1 = " << toString(other.X) << "*" << toString(Z) << " = " << toString(V1) << mp << std::endl;
            std::cout << "  V2 = X1*Z2 = " << toString(X) << "*" << toString(other.Z) << " = " << toString(V2) << mp << std::endl;
            std::cout << "  U = U1-U2 = " << toString(U1) << "-" << toString(U2) << " = " << toString(U) << mp << std::endl;
            std::cout << "  V = V1-V2 = " << toString(V1) << "-" << toString(V2) << " = " << toString(V) << mp << std::endl;
            std::cout << "  W = Z1*Z2 = " << toString(Z) << "*" << toString(other.Z) << " = " << toString(W) << mp << std::endl;
            std::cout << "  A = U^2*W - V^3 - 2*V^2*V2 = " << toString(A) << mp << std::endl;
            std::cout << "  X3 = V*A = " << toString(V) << "*" << toString(A) << " = " << toString(X3) << mp << std::endl;
            std::cout << "  Y3 = U*(V^2*V2-A) - V^3*U2 = " << toString(Y3) << mp << std::endl;
            std::cout << "  Z3 = V^3*W = " << toString(Vcb) << "*" << toString(W) << " = " << toString(Z3) << mp << std::endl;
        }

        return EllipticCurvePoint(X3, Y3, Z3, curve);
    }

    // --- Scalar Multiplication ---

    // Double-and-Add algorithm: kP = P + P + ... + P (k times)
    // Iterates over bits of k from LSB to MSB.
    // Running time depends on the number of set bits (NOT constant-time).
    EllipticCurvePoint scalarMul(const T& k) const {
        if (k == 0 || isInfinity())
            return infinity(curve);

        assertOnCurve("scalarMul");

        EllipticCurvePoint res = infinity(curve);   // accumulator
        EllipticCurvePoint temp = *this;            // current power: P, 2P, 4P, ...

        std::vector<int> bits = getBits(k);

        if (curve->verbose) {
            std::cout << "[scalarMul DoubleAndAdd] k = " << toString(k)
                      << " (" << bits.size() << " bits)" << std::endl;
        }

        for (size_t i = 0; i < bits.size(); i++) {
            if (bits[i] == 1) {
                if (curve->verbose)
                    std::cout << "  bit[" << i << "] = 1: res = res + temp" << std::endl;
                res = res.pointAdd(temp);
            } else {
                if (curve->verbose)
                    std::cout << "  bit[" << i << "] = 0: skip" << std::endl;
            }
            temp = temp.pointDouble();
        }

        return res;
    }

    // Montgomery ladder algorithm: kP
    // Constant-time execution — always performs one addition and one doubling per bit.
    // Important for cryptography (protection against side-channel attacks).
    EllipticCurvePoint scalarMulMontgomery(const T& k) const {
        if (k == 0 || isInfinity())
            return infinity(curve);

        assertOnCurve("scalarMulMontgomery");

        EllipticCurvePoint R0 = infinity(curve);
        EllipticCurvePoint R1 = *this;

        std::vector<int> bits = getBits(k);

        if (curve->verbose) {
            std::cout << "[scalarMul Montgomery] k = " << toString(k)
                      << " (" << bits.size() << " bits)" << std::endl;
        }

        // Reverse iteration: from MSB to LSB
        for (int i = static_cast<int>(bits.size()) - 1; i >= 0; i--) {
            if (bits[i] == 0) {
                if (curve->verbose)
                    std::cout << "  bit[" << i << "] = 0: R1 = R0+R1, R0 = 2*R0" << std::endl;
                R1 = R0.pointAdd(R1);
                R0 = R0.pointDouble();
            } else {
                if (curve->verbose)
                    std::cout << "  bit[" << i << "] = 1: R0 = R0+R1, R1 = 2*R1" << std::endl;
                R0 = R0.pointAdd(R1);
                R1 = R1.pointDouble();
            }
        }

        return R0;
    }

private:
    // Binary representation of a number (LSB first).
    // Works for both long long and mpz_class (both support % and /).
    static std::vector<int> getBits(const T& k) {
        std::vector<int> bits;
        T val = k;
        while (val > 0) {
            T remainder = val % 2;
            bits.push_back(remainder == 0 ? 0 : 1);
            val /= 2;
        }
        return bits;  // bits[0] = LSB, bits[size-1] = MSB
    }

public:

    // --- Point Order ---

    // Finds the order of point P — the smallest k > 0 such that kP = O_E.
    // By Lagrange's theorem, the point order divides the curve order n.
    // Therefore we enumerate divisors of n from smallest to largest.
    T pointOrder() const {
        if (isInfinity()) return T(1);

        const T& n = curve->n;
        if (n == 0) throw std::runtime_error("pointOrder: curve order n is not set");

        // Collect divisors of n (set = automatically sorted)
        std::set<T> divisors;
        for (T d = 1; d * d <= n; d += 1) {
            if (n % d == 0) {
                divisors.insert(d);
                divisors.insert(n / d);
            }
        }

        // Check from smallest
        for (const T& d : divisors) {
            auto res = scalarMul(d);
            if (res.isInfinity()) {
                if (curve->verbose)
                    std::cout << "[pointOrder] ord(P) = " << toString(d) << " (first divisor of n=" << toString(n) << " where d*P = O_E)" << std::endl;
                return d;
            }
        }

        throw std::runtime_error("pointOrder: failed to find order (should not happen)");
    }

    // --- Output ---

    // Helper: T -> string (for MSVC compatibility with mpz_class operator<<)
    static std::string toString(const T& val) {
        std::ostringstream oss;
        oss << val;
        return oss.str();
    }

    // Print the point. Verbose does NOT affect this — output is always compact.
    void print() const {
        if (isInfinity()) {
            std::cout << "O_E (0 : 1 : 0)" << std::endl;
        } else {
            // Compute affine coordinates silently (without verbose)
            const T& p = curve->p;
            T zInv = modInverse(Z, p);
            T ax = mod(X * zInv, p);
            T ay = mod(Y * zInv, p);
            std::cout << "(" << toString(ax) << ", " << toString(ay) << ")  "
                      << "[projective: (" << toString(X) << " : " << toString(Y) << " : " << toString(Z) << ")]"
                      << std::endl;
        }
    }

    // --- Comparison ---
    // Two projective points are equal if (X1*Z2 == X2*Z1) and (Y1*Z2 == Y2*Z1)
    bool equals(const EllipticCurvePoint& other) const {
        if (isInfinity() && other.isInfinity()) return true;
        if (isInfinity() || other.isInfinity()) return false;

        const T& p = curve->p;
        return mod(X * other.Z, p) == mod(other.X * Z, p) &&
               mod(Y * other.Z, p) == mod(other.Y * Z, p);
    }
};
