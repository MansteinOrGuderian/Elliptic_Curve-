#pragma once

// Еліптична крива у формі Вейєрштраса: y^2 = x^3 + ax + b (mod p)
// T = mpz_class для Baby-JubJub, T = long long для спрощених випадків

template <typename T>
class EllipticCurve {
public:
    T a;        // коефіцієнт a
    T b;        // коефіцієнт b
    T p;        // модуль поля (просте число)
    T n;        // порядок кривої (кількість точок, включно з O_E)
    bool verbose;  // покроковий вивід обчислень (для long long)

    // Конструктор з усіма параметрами
    EllipticCurve(const T& a, const T& b, const T& p, const T& n, bool verbose = false)
        : a(a), b(b), p(p), n(n), verbose(verbose) {}

    // Допоміжна функція: T -> string
    static std::string toString(const T& val) {
        std::ostringstream oss;
        oss << val;
        return oss.str();
    }

    // Перевірка що крива невироджена: 4a^3 + 27b^2 != 0 (mod p)
    bool isNonSingular() const {
        T discriminant = mod(4 * a * a * a + 27 * b * b, p);
        if (verbose) {
            std::cout << "[isNonSingular] 4*" << toString(a) << "^3 + 27*" << toString(b) << "^2 = "
                      << toString(discriminant) << " (mod " << toString(p) << ") "
                      << (discriminant != 0 ? "!= 0 -> curve is non-singular" : "== 0 -> curve is SINGULAR")
                      << std::endl;
        }
        return discriminant != 0;
    }

    // Виведення інформації про криву
    void print() const {
        std::cout << "Elliptic Curve: y^2 = x^3";
        if (a != 0) std::cout << " + " << toString(a) << "*x";
        if (b != 0) std::cout << " + " << toString(b);
        std::cout << "  (mod " << toString(p) << ")" << std::endl;
        std::cout << "Order n = " << toString(n) << std::endl;
    }
};
