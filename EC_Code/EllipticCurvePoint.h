#pragma once

// Точка еліптичної кривої у проективних координатах (X, Y, Z)
// Афінна точка (x, y) відповідає проективній (x, y, 1)
// Точка на нескінченності O_E = (0, 1, 0)

template <typename T>
class EllipticCurvePoint {
public:
    T X, Y, Z;
    const EllipticCurve<T>* curve;  // вказівник на криву, якій належить точка

    // --- Конструктори ---

    // Конструктор проективної точки
    EllipticCurvePoint(const T& X, const T& Y, const T& Z, const EllipticCurve<T>* curve)
        : X(X), Y(Y), Z(Z), curve(curve) {}

    // Створення точки на нескінченності
    static EllipticCurvePoint infinity(const EllipticCurve<T>* curve) {
        return EllipticCurvePoint(T(0), T(1), T(0), curve);
    }

    // Створення точки з афінних координат (x, y) -> (x, y, 1)
    static EllipticCurvePoint fromAffine(const T& x, const T& y, const EllipticCurve<T>* curve) {
        return EllipticCurvePoint(x, y, T(1), curve);
    }

    // --- Перевірки ---

    bool isInfinity() const {
        return Z == 0;
    }

    // Перевірка належності точки до кривої: Y^2*Z = X^3 + a*X*Z^2 + b*Z^3 (mod p)
    bool isOnCurve() const {
        if (isInfinity()) return true;

        const T& p = curve->p;
        const T& a = curve->a;
        const T& b = curve->b;

        // mod після кожного множення, щоб уникнути переповнення long long
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

    // --- Конвертація ---

    // Проективні -> афінні: (X, Y, Z) -> (X/Z, Y/Z)
    // Повертає пару (x, y). Для O_E кидає виключення.
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

    // --- Арифметика ---

    // Подвоєння точки у проективних координатах
    // За псевдокодом: W = a*Z^2 + 3*X^2, S = Y*Z, B = X*Y*S,
    //   H = W^2 - 8*B, X' = 2*H*S, Y' = W*(4*B - H) - 8*Y^2*S^2, Z' = 8*S^3
    EllipticCurvePoint pointDouble() const {
        // 2 * O_E = O_E
        if (isInfinity()) {
            if (curve->verbose)
                std::cout << "[pointDouble] O_E doubled -> O_E" << std::endl;
            return infinity(curve);
        }

        const T& p = curve->p;
        const T& a = curve->a;

        // Точка порядку 2: Y == 0 => 2P = O_E
        if (mod(Y, p) == 0) {
            if (curve->verbose)
                std::cout << "[pointDouble] Y == 0, point of order 2 -> O_E" << std::endl;
            return infinity(curve);
        }

        // Подвоєння у проективних координатах.
        // Формули з псевдокоду: W, S, B, H -> X', Y', Z'.
        // mod() після кожного множення щоб уникнути переповнення long long.
        T W  = mod(mod(a * mod(Z * Z, p), p) + mod(3 * mod(X * X, p), p), p);  // W = a*Z^2 + 3*X^2
        T S  = mod(Y * Z, p);                                                   // S = Y*Z
        T B  = mod(mod(X * Y, p) * S, p);                                       // B = X*Y*S
        T H  = mod(mod(W * W, p) - mod(8 * B, p), p);                           // H = W^2 - 8*B
        T Xr = mod(mod(2 * H, p) * S, p);                                       // X' = 2*H*S
        T Yr = mod(mod(W * mod(4 * B - H, p), p)                                // Y' = W*(4B-H)
                  - mod(mod(8 * mod(Y * Y, p), p) * mod(S * S, p), p), p);      //     - 8*Y^2*S^2
        T Zr = mod(8 * mod(mod(S * S, p) * S, p), p);                           // Z' = 8*S^3

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

    // --- Виведення ---

    // Допоміжна функція: T -> string (для сумісності з MSVC operator<< для mpz_class)
    static std::string toString(const T& val) {
        std::ostringstream oss;
        oss << val;
        return oss.str();
    }

    // Виводить точку. Verbose тут НЕ впливає — вивід завжди компактний.
    void print() const {
        if (isInfinity()) {
            std::cout << "O_E (0 : 1 : 0)" << std::endl;
        } else {
            // Обчислюємо афінні координати тихо (без verbose)
            const T& p = curve->p;
            T zInv = modInverse(Z, p);
            T ax = mod(X * zInv, p);
            T ay = mod(Y * zInv, p);
            std::cout << "(" << toString(ax) << ", " << toString(ay) << ")  "
                      << "[projective: (" << toString(X) << " : " << toString(Y) << " : " << toString(Z) << ")]"
                      << std::endl;
        }
    }

    // --- Порівняння ---
    // Дві проективні точки рівні якщо (X1*Z2 == X2*Z1) та (Y1*Z2 == Y2*Z1)
    bool equals(const EllipticCurvePoint& other) const {
        if (isInfinity() && other.isInfinity()) return true;
        if (isInfinity() || other.isInfinity()) return false;

        const T& p = curve->p;
        return mod(X * other.Z, p) == mod(other.X * Z, p) &&
               mod(Y * other.Z, p) == mod(other.Y * Z, p);
    }
};
