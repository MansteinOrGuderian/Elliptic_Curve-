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

    // Перевірка що крива невироджена: Δ = -16(4a^3 + 27b^2) != 0 (mod p)
    bool isNonSingular() const {
        T disc_core = mod(4 * a * a * a + 27 * b * b, p);
        if (verbose) {
            std::cout << "[isNonSingular] delta = -16*(4a^3 + 27b^2): 4*" << toString(a) << "^3 + 27*" << toString(b) << "^2 = "
                      << toString(disc_core) << " (mod " << toString(p) << ") "
                      << (disc_core != 0 ? "!= 0 -> non-singular" : "== 0 -> SINGULAR!")
                      << std::endl;
        }
        return disc_core != 0;
    }

    // Виведення рівняння кривої
    void print() const {
        std::cout << "Elliptic Curve: y^2 = x^3";
        if (a != 0) std::cout << " + " << toString(a) << "*x";
        if (b != 0) std::cout << " + " << toString(b);
        std::cout << "  (mod " << toString(p) << ")" << std::endl;
    }

    // Знаходження УСІХ точок кривої перебором (тільки для малих p!)
    // Зберігає порядок в n. Повертає вектор афінних точок (без O_E).
    // sqrtExampleX: для якого x показати покрокове обчислення кореня (-1 = не показувати).
    std::vector<std::pair<long long, long long>> findAllPoints(bool showTable = true, long long sqrtExampleX = 0) {
        long long pp = static_cast<long long>(p);

        if (verbose) {
            if (pp % 4 == 3)
                std::cout << "p = " << pp << " = 4k+3 (k=" << (pp - 3) / 4 << "), sqrt via y = a^(k+1)" << std::endl;
            else if (pp % 8 == 5)
                std::cout << "p = " << pp << " = 8k+5 (k=" << (pp - 5) / 8 << ")" << std::endl;
            else
                std::cout << "p = " << pp << " = 8k+1, using Tonelli-Shanks" << std::endl;
        }

        std::vector<long long> rhs(pp);
        std::vector<long long> y1arr(pp, -1);
        std::vector<std::pair<long long, long long>> points;

        long long count = 1;  // O_E

        for (long long x = 0; x < pp; x++) {
            rhs[x] = mod(modPow(x, 3LL, pp) + mod(static_cast<long long>(a) * x, pp) + static_cast<long long>(b), pp);

            if (rhs[x] == 0) {
                y1arr[x] = 0;
                points.push_back({x, 0});
                count++;
            } else {
                long long y = modSqrt(rhs[x], pp, /*verbose=*/false);
                if (y >= 0) {
                    y1arr[x] = y;
                    points.push_back({x, y});
                    points.push_back({x, mod(-y, pp)});
                    count += 2;
                }
            }
        }

        n = static_cast<T>(count);

        // Покрокове обчислення кореня для обраного x
        if (verbose && sqrtExampleX >= 0 && sqrtExampleX < pp) {
            std::cout << "\nSqrt example for x=" << sqrtExampleX
                      << ": RHS = " << rhs[sqrtExampleX] << " (mod " << pp << ")" << std::endl;
            if (rhs[sqrtExampleX] == 0) {
                std::cout << "  RHS = 0 -> y = 0 (single point)" << std::endl;
            } else {
                modSqrt(rhs[sqrtExampleX], pp, /*verbose=*/true);
                if (y1arr[sqrtExampleX] > 0)
                    std::cout << "  => y1 = " << y1arr[sqrtExampleX] << ", y2 = " << mod(-y1arr[sqrtExampleX], pp) << std::endl;
            }
            std::cout << std::endl;
        }

        // Таблиця (афінні координати)
        if (showTable) {
            const int MAX_COLS = 15;     // макс. колонок в одному рядку таблиці
            const int LABEL_W  = 7;      // ширина колонки підписів (x, RHS, y1, y2)

            std::cout << "--- All " << count << " affine points of y^2 = x^3 + "
                      << toString(a) << "*x + " << toString(b) << " (mod " << toString(p) << ") ---\n";
            std::cout << "(+ O_E, which has no affine representation)\n\n";

            for (long long start = 0; start < pp; start += MAX_COLS) {
                long long end = std::min(start + MAX_COLS, pp);
                int ncols = static_cast<int>(end - start);
                int sepLen = LABEL_W + 1 + 8 * ncols;

                auto sep = [&]() {
                    std::cout << "|";
                    for (int i = 0; i < sepLen - 1; i++) std::cout << "-";
                    std::cout << "|\n";
                };

                sep();
                std::cout << "|" << std::setw(LABEL_W) << "x" << "|";
                for (long long x = start; x < end; x++)
                    std::cout << std::setw(6) << x << " |";
                std::cout << "\n";
                sep();

                std::cout << "|" << std::setw(LABEL_W) << "RHS" << "|";
                for (long long x = start; x < end; x++)
                    std::cout << std::setw(6) << rhs[x] << " |";
                std::cout << "\n";
                sep();

                std::cout << "|" << std::setw(LABEL_W) << "y1" << "|";
                for (long long x = start; x < end; x++) {
                    if (y1arr[x] >= 0)
                        std::cout << std::setw(6) << y1arr[x] << " |";
                    else
                        std::cout << std::setw(6) << "-" << " |";
                }
                std::cout << "\n";
                sep();

                std::cout << "|" << std::setw(LABEL_W) << "y2" << "|";
                for (long long x = start; x < end; x++) {
                    if (y1arr[x] > 0)
                        std::cout << std::setw(6) << mod(-y1arr[x], pp) << " |";
                    else
                        std::cout << std::setw(6) << "-" << " |";
                }
                std::cout << "\n";
                sep();
                std::cout << "\n";
            }
        }

        return points;
    }
};
