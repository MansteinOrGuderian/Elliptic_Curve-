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
    // Перевіряє належність точки до кривої; кидає виключення якщо ні.
    static EllipticCurvePoint fromAffine(const T& x, const T& y, const EllipticCurve<T>* curve) {
        EllipticCurvePoint pt(x, y, T(1), curve);
        if (!pt.isOnCurve())
            throw std::runtime_error("fromAffine: point (" + toString(x) + ", " + toString(y) + ") is NOT on the curve");
        return pt;
    }

    // Створення точки БЕЗ перевірки (для внутрішнього використання —
    // результати арифметики гарантовано на кривій, якщо вхід був коректний)
    static EllipticCurvePoint fromAffineUnchecked(const T& x, const T& y, const EllipticCurve<T>* curve) {
        return EllipticCurvePoint(x, y, T(1), curve);
    }

    // --- Перевірки ---

    bool isInfinity() const {
        return Z == 0;
    }

    // Перевірка з throw — для використання перед арифметичними операціями.
    // Тиха: не друкує verbose, бо це внутрішня перевірка.
    void assertOnCurve(const std::string& context) const {
        if (isInfinity()) return;  // O_E завжди на кривій

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

        assertOnCurve("pointDouble");

        const T& p = curve->p;
        const T& a = curve->a;

        // Точка порядку 2: Y == 0 => 2P = O_E
        if (mod(Y, p) == 0) {
            if (curve->verbose)
                std::cout << "[pointDouble] P = (" << toString(X) << ", " << toString(Y) << ", " << toString(Z)
                          << "), Y == 0 -> point of order 2 -> O_E" << std::endl;
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

    // Додавання двох точок у проективних координатах
    // За псевдокодом: U1, U2, V1, V2 -> U, V, W, A -> X3, Y3, Z3
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

        // U1 = Y2*Z1,  U2 = Y1*Z2  (порівняння Y-координат у спільному масштабі)
        T U1 = mod(other.Y * Z, p);
        T U2 = mod(Y * other.Z, p);
        // V1 = X2*Z1,  V2 = X1*Z2  (порівняння X-координат у спільному масштабі)
        T V1 = mod(other.X * Z, p);
        T V2 = mod(X * other.Z, p);

        if (V1 == V2) {
            // Однакова X-координата
            if (U1 != U2) {
                // P і Q — взаємно обернені: P + (-P) = O_E
                if (curve->verbose)
                    std::cout << "[pointAdd] P = (" << toString(X) << ", " << toString(Y) << ", " << toString(Z)
                              << "), Q = (" << toString(other.X) << ", " << toString(other.Y) << ", " << toString(other.Z)
                              << ") -> inverse points -> O_E" << std::endl;
                return infinity(curve);
            } else {
                // P == Q: переходимо до подвоєння
                if (curve->verbose)
                    std::cout << "[pointAdd] P == Q -> calling pointDouble" << std::endl;
                return pointDouble();
            }
        }

        // Загальний випадок: P != Q, P != -Q
        T U  = mod(U1 - U2, p);                          // U = U1 - U2
        T V  = mod(V1 - V2, p);                          // V = V1 - V2
        T W  = mod(Z * other.Z, p);                      // W = Z1*Z2
        T Vsq = mod(V * V, p);                           // V^2
        T Vcb = mod(Vsq * V, p);                         // V^3
        T Usq = mod(U * U, p);                           // U^2
        T VsqV2 = mod(Vsq * V2, p);                     // V^2 * V2
        T A  = mod(mod(Usq * W, p) - Vcb - mod(2 * VsqV2, p), p);  // A = U^2*W - V^3 - 2*V^2*V2
        T X3 = mod(V * A, p);                            // X3 = V*A
        T Y3 = mod(mod(U * mod(VsqV2 - A, p), p) - mod(Vcb * U2, p), p);  // Y3 = U*(V^2*V2 - A) - V^3*U2
        T Z3 = mod(Vcb * W, p);                          // Z3 = V^3*W

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

    // --- Скалярний добуток ---

    // Алгоритм DoubleAndAdd: kP = P + P + ... + P (k разів)
    // Ітерує по бітах k від молодшого до старшого.
    // Час виконання залежить від кількості одиничних бітів (НЕ константний).
    EllipticCurvePoint scalarMul(const T& k) const {
        if (k == 0 || isInfinity())
            return infinity(curve);

        assertOnCurve("scalarMul");

        EllipticCurvePoint res = infinity(curve);   // акумулятор
        EllipticCurvePoint temp = *this;            // поточна степінь: P, 2P, 4P, ...

        // Отримуємо бітове представлення k
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

    // Алгоритм сходів Монтгомері: kP
    // Константний час виконання — завжди робить одне додавання і одне подвоєння на біт.
    // Це важливо для криптографії (захист від side-channel атак).
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

        // Зворотня ітерація: від старшого біта до молодшого
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
    // Бітове представлення числа (від молодшого до старшого біту).
    // Працює і для long long, і для mpz_class (обидва підтримують % і /).
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

    // --- Порядок точки ---

    // Знаходить порядок точки P — найменше k > 0 таке що kP = O_E.
    // За теоремою Лагранжа, порядок точки ділить порядок кривої n.
    // Тому перебираємо дільники n від найменшого.
    T pointOrder() const {
        if (isInfinity()) return T(1);

        const T& n = curve->n;
        if (n == 0) throw std::runtime_error("pointOrder: curve order n is not set");

        // Збираємо дільники n
        std::vector<T> divisors;
        for (T d = 1; d * d <= n; d += 1) {
            if (n % d == 0) {
                divisors.push_back(d);
                if (d != n / d) divisors.push_back(n / d);
            }
        }
        std::sort(divisors.begin(), divisors.end());

        // Перевіряємо від найменшого
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
