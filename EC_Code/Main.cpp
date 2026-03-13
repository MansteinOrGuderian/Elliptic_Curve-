#include "Header.h"

int main() {
    std::cout << "========== Simplified case ==========\n\n";
    {
        // Edit a, b, p кривої y^2 = x^3 + a*x + b (mod p)
        long long sa = 5, sb = 7, sp = 11;

        EllipticCurve<long long> curve(sa, sb, sp, 0LL, /*verbose=*/true);
        curve.print();
        curve.isNonSingular();
        std::cout << std::endl;

        // Знаходимо всі точки (заповнює n); sqrtExampleX=2 — показати обчислення для x=2
        auto allPoints = curve.findAllPoints(/*showTable=*/true, /*sqrtExampleX=*/2);
        std::cout << "Order of curve: n = " << curve.n << std::endl;
        std::cout << std::endl;

        // Беремо першу знайдену точку з y != 0 як P
        // (точки з y=0 мають порядок 2 — для демо цікавіші "звичайні" точки)
        long long px = 0, py = 0;
        for (auto& [x, y] : allPoints) {
            if (y != 0) { px = x; py = y; break; }
        }
        auto P = EllipticCurvePoint<long long>::fromAffine(px, py, &curve);
        std::cout << "P = ";
        P.print();
        std::cout << std::endl;

        // Спроба створити точку з неправильним y — fromAffine кине виключення
        long long fakeY = mod(py + 1, static_cast<long long>(curve.p));
        if (fakeY == mod(-py, static_cast<long long>(curve.p))) fakeY = mod(py + 2, static_cast<long long>(curve.p));  // щоб не потрапити на -P
        try {
            auto bad = EllipticCurvePoint<long long>::fromAffine(px, fakeY, &curve);
        } catch (const std::runtime_error& e) {
            std::cout << e.what() << std::endl;
        }
        std::cout << std::endl;

        // --- PointDouble тест ---
        std::cout << "--- PointDouble test ---" << std::endl;
        auto twoP = P.pointDouble();
        std::cout << "2P = ";
        twoP.print();
        std::cout << "2P on curve: " << (twoP.isOnCurve() ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // 2*(7,0): Y=0 => точка порядку 2, має дати O_E
        auto fourP = twoP.pointDouble();
        std::cout << "4P = 2*(2P) = ";
        fourP.print();
        std::cout << std::endl;

        // --- PointAdd тести ---
        std::cout << "--- PointAdd tests ---" << std::endl;

        // Q — друга знайдена точка з y != 0 і x != px
        long long qx = 0, qy = 0;
        for (auto& [x, y] : allPoints) {
            if (y != 0 && x != px) { qx = x; qy = y; break; }
        }
        auto Q = EllipticCurvePoint<long long>::fromAffine(qx, qy, &curve);
        std::cout << "Q = ";
        Q.print();

        // -P = (px, -py mod p)
        long long negPy = mod(-py, static_cast<long long>(curve.p));
        auto negP = EllipticCurvePoint<long long>::fromAffine(px, negPy, &curve);
        std::cout << "-P = ";
        negP.print();
        std::cout << std::endl;

        // P + Q (різні точки)
        auto PpQ = P.pointAdd(Q);
        std::cout << "P + Q = ";
        PpQ.print();
        std::cout << "On curve: " << (PpQ.isOnCurve() ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // P + P має дорівнювати 2P (перевірка що pointAdd делегує в pointDouble)
        auto PpP = P.pointAdd(P);
        std::cout << "P + P = ";
        PpP.print();
        std::cout << "Equals 2P: " << (PpP.equals(twoP) ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // P + (-P) = O_E
        auto PpNegP = P.pointAdd(negP);
        std::cout << "P + (-P) = ";
        PpNegP.print();
        std::cout << std::endl;

        // P + O_E = P
        auto O_E = EllipticCurvePoint<long long>::infinity(&curve);
        auto PpO = P.pointAdd(O_E);
        std::cout << "P + O_E = ";
        PpO.print();
        std::cout << "Equals P: " << (PpO.equals(P) ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // --- ScalarMul тести ---
        // Вимикаємо verbose — покрокову арифметику вже бачили вище
        curve.verbose = false;
        std::cout << "--- ScalarMul tests (verbose off) ---" << std::endl;

        // 3P двома алгоритмами
        auto threeP_daa = P.scalarMul(3LL);
        std::cout << "3P (DoubleAndAdd) = ";
        threeP_daa.print();

        auto threeP_mont = P.scalarMulMontgomery(3LL);
        std::cout << "3P (Montgomery)   = ";
        threeP_mont.print();

        std::cout << "Both equal: " << (threeP_daa.equals(threeP_mont) ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // Перевірка: 3P == P + 2P
        auto manual3P = P.pointAdd(twoP);
        std::cout << "P + 2P = ";
        manual3P.print();
        std::cout << "3P == P + 2P: " << (threeP_daa.equals(manual3P) ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // Крайні випадки
        auto zeroP = P.scalarMul(0LL);
        std::cout << "0*P = ";
        zeroP.print();

        auto oneP = P.scalarMul(1LL);
        std::cout << "1*P = ";
        oneP.print();
        std::cout << "Equals P: " << (oneP.equals(P) ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // Порядок точки P
        long long ordP = P.pointOrder();
        std::cout << "ord(P) = " << ordP << "  (divides n=" << curve.n << ")" << std::endl;
        std::cout << std::endl;

        // Перевірка пункту 2 завдання: n*P = O_E (n — порядок кривої)
        auto nP = P.scalarMul(curve.n);
        std::cout << "Curve order verification: n*P (n=" << curve.n << ") = ";
        nP.print();
    }

    std::cout << "\n========== Baby-JubJub ==========\n\n";
    {
        mpz_class p_bjj("21888242871839275222246405745257275088548364400416034343698204186575808495617");
        std::cout << "Field: F_p,  p = " << p_bjj << std::endl;
        std::cout << std::endl;

        // 1) Twisted Edwards form: a_e*x^2 + y^2 = 1 + d*x^2*y^2
        mpz_class a_edwards(168700);
        mpz_class d_edwards("9706598848417545097372247223557719406784115219466060233080913168975159366771");
        std::cout << "Twisted Edwards: " << a_edwards << "*x^2 + y^2 = 1 + d*x^2*y^2" << std::endl;
        std::cout << "  d = " << d_edwards << std::endl;
        std::cout << std::endl;

        // 2) Montgomery form: v^2 = u^3 + A*u^2 + u
        mpz_class A_mont(168698);
        std::cout << "Montgomery: v^2 = u^3 + " << A_mont << "*u^2 + u" << std::endl;
        std::cout << std::endl;

        // 3) Weierstrass form: y^2 = x^3 + a*x + b
        // Конвертація: a = (3 - A^2)/3,  b = (2A^3 - 9A)/27
        mpz_class inv3  = modInverse(mpz_class(3), p_bjj);
        mpz_class inv27 = modInverse(mpz_class(27), p_bjj);
        mpz_class A2 = mod(A_mont * A_mont, p_bjj);
        mpz_class A3 = mod(A2 * A_mont, p_bjj);
        mpz_class a_bjj = mod((3 - A2) * inv3, p_bjj);
        mpz_class b_bjj = mod((2 * A3 - 9 * A_mont) * inv27, p_bjj);

        // Порядок: n = 8*r (з документації)
        mpz_class r_bjj("2736030358979909402780800718157159386076813972158567259200215660948447373041");
        mpz_class n_bjj = 8 * r_bjj;

        EllipticCurve<mpz_class> bjj(a_bjj, b_bjj, p_bjj, n_bjj, /*verbose=*/false);
        std::cout << "Weierstrass: ";
        bjj.print();
        std::cout << "Non-singular: " << (bjj.isNonSingular() ? "YES" : "NO") << std::endl;
        std::cout << "Order n = 8 * r = " << n_bjj << std::endl;
        std::cout << "  r = " << r_bjj << " (prime, subgroup order)" << std::endl;
        std::cout << std::endl;

        // Генератор: Edwards -> Montgomery -> Weierstrass
        // Edwards generator (з документації iden3):
        mpz_class ex("5299619240641551281634865583518297030282874472190772894086521144482721001553");
        mpz_class ey("16950150798460657717958625567821834550301663161624707787222815936182638968203");

        // Edwards -> Montgomery: u = (1+y)/(1-y),  v = u/x
        mpz_class u = mod((1 + ey) * modInverse(mod(1 - ey, p_bjj), p_bjj), p_bjj);
        mpz_class v = mod(u * modInverse(ex, p_bjj), p_bjj);

        // Montgomery -> Weierstrass: x_w = u + A/3,  y_w = v
        mpz_class gx = mod(u + A_mont * inv3, p_bjj);
        mpz_class gy = v;

        auto G = EllipticCurvePoint<mpz_class>::fromAffine(gx, gy, &bjj);
        std::cout << "Generator G: ";
        G.print();

        // PointDouble / PointAdd тести
        auto twoG = G.pointDouble();
        std::cout << "2G = ";
        twoG.print();

        auto GpG = G.pointAdd(G);
        std::cout << "G + G equals 2G: " << (GpG.equals(twoG) ? "YES" : "NO") << std::endl;

        auto threeG = G.pointAdd(twoG);
        std::cout << "3G = ";
        threeG.print();
        std::cout << "3G on curve: " << (threeG.isOnCurve() ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // ScalarMul тести
        std::cout << "--- ScalarMul tests (BJJ) ---" << std::endl;

        auto fiveG_daa  = G.scalarMul(mpz_class(5));
        auto fiveG_mont = G.scalarMulMontgomery(mpz_class(5));
        std::cout << "5G (DoubleAndAdd) = ";
        fiveG_daa.print();
        std::cout << "5G (Montgomery)   = ";
        fiveG_mont.print();
        std::cout << "Both equal: " << (fiveG_daa.equals(fiveG_mont) ? "YES" : "NO") << std::endl;
        std::cout << std::endl;

        // Ключова перевірка (пункт 2 завдання): n*G = O_E
        std::cout << "Verifying n*G = O_E..." << std::endl;

        auto start_daa = std::chrono::high_resolution_clock::now();
        auto nG_daa = G.scalarMul(n_bjj);
        auto end_daa = std::chrono::high_resolution_clock::now();
        double ms_daa = std::chrono::duration<double, std::milli>(end_daa - start_daa).count();
        std::cout << "n*G (DoubleAndAdd):  " << (nG_daa.isInfinity() ? "O_E" : "NOT O_E")
                  << "  [" << ms_daa << " ms]" << std::endl;

        auto start_mont = std::chrono::high_resolution_clock::now();
        auto nG_mont = G.scalarMulMontgomery(n_bjj);
        auto end_mont = std::chrono::high_resolution_clock::now();
        double ms_mont = std::chrono::duration<double, std::milli>(end_mont - start_mont).count();
        std::cout << "n*G (Montgomery):    " << (nG_mont.isInfinity() ? "O_E" : "NOT O_E")
                  << "  [" << ms_mont << " ms]" << std::endl;
        std::cout << std::endl;
    }

    std::cout << "\n========== Modular arithmetic tests ==========\n\n";
    {
        std::cout << "mod(-3, 11) = " << mod(-3LL, 11LL) << " (expected 8)" << std::endl;
        std::cout << "modInverse(3, 11) = " << modInverse(3LL, 11LL) << " (expected 4)" << std::endl;
        std::cout << "modPow(2, 10, 11) = " << modPow(2LL, 10LL, 11LL) << " (expected 1)" << std::endl;
    }

    return 0;
}
