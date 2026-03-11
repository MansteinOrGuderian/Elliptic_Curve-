#include "Header.h"

int main() {
    std::cout << "========== Simplified case: y^2 = x^3 + 5x + 7 (mod 11) ==========\n\n";
    {
        EllipticCurve<long long> curve(5LL, 7LL, 11LL, 0LL, /*verbose=*/true);
        curve.print();
        std::cout << std::endl;

        // isNonSingular: викликаємо окремо, щоб verbose вивід не вклинювався
        curve.isNonSingular();
        std::cout << std::endl;

        // Точка на нескінченності
        auto O = EllipticCurvePoint<long long>::infinity(&curve);
        std::cout << "O_E: ";
        O.print();
        std::cout << std::endl;

        // Точка (4, 5) -- на кривій:
        // x=4: x^3+5x+7 = 64+20+7 = 91 = 3 (mod 11), y=5: y^2 = 25 = 3 (mod 11) -> OK
        auto P = EllipticCurvePoint<long long>::fromAffine(4LL, 5LL, &curve);
        std::cout << "Point P: ";
        P.print();
        P.isOnCurve();
        std::cout << std::endl;

        // Точка (1, 1) -- НЕ на кривій
        auto fake = EllipticCurvePoint<long long>::fromAffine(1LL, 1LL, &curve);
        std::cout << "Fake point: ";
        fake.print();
        fake.isOnCurve();
        std::cout << std::endl;

        // --- PointDouble тест ---
        // 2P для P=(4,5): очікуємо (7, 0) за ручним обрахунком
        // λ = (3*16+5)/(2*5) = 53/10 = 9*10^{-1} = 9*10 = 90 = 2 (mod 11)
        // x3 = 4 - 8 = -4 = 7 (mod 11), y3 = 2*(4-7) - 5 = -11 = 0 (mod 11)
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
    }

    std::cout << "\n========== Baby-JubJub ==========\n\n";
    {
        mpz_class p_bjj("21888242871839275222246405745257275088548364400416034343698204186575808495617");
        mpz_class a_bjj("7296080957279758407415468581752425029516121466805344781232734728849116493472");
        mpz_class b_bjj("16213513238399463127589930181672055621146936592900766180517188641980520820846");
        mpz_class n_bjj("21888242871839275222246405745257275088614511777268538073601725287587578984328");

        EllipticCurve<mpz_class> bjj(a_bjj, b_bjj, p_bjj, n_bjj, /*verbose=*/false);
        bjj.print();
        std::cout << "Non-singular: " << (bjj.isNonSingular() ? "YES" : "NO") << std::endl;

        auto O = EllipticCurvePoint<mpz_class>::infinity(&bjj);
        std::cout << "O_E: ";
        O.print();

        // Генератор Baby-JubJub: конвертуємо з Едвардса -> Монтгомері -> Вейєрштрас
        // Twisted Edwards: 168700*x^2 + y^2 = 1 + 168696*x^2*y^2
        // Edwards generator:
        mpz_class ex("5299619240641551281634865583518297030282874472190772894086521144482721001553");
        mpz_class ey("16950150798460657717958625567821834550301663161624707787222815936182638968203");

        // Edwards -> Montgomery: u = (1+y)/(1-y),  v = u/x
        mpz_class u = mod((1 + ey) * modInverse(mod(1 - ey, p_bjj), p_bjj), p_bjj);
        mpz_class v = mod(u * modInverse(ex, p_bjj), p_bjj);

        // Montgomery -> Weierstrass: x_w = u + A/3,  y_w = v
        mpz_class A_mont(168698);
        mpz_class inv3 = modInverse(mpz_class(3), p_bjj);
        mpz_class gx = mod(u + A_mont * inv3, p_bjj);
        mpz_class gy = v;

        auto G = EllipticCurvePoint<mpz_class>::fromAffine(gx, gy, &bjj);
        std::cout << "Generator G: ";
        G.print();
        std::cout << "G is on curve: " << (G.isOnCurve() ? "YES" : "NO") << std::endl;

        // --- PointDouble тест для BJJ ---
        auto twoG = G.pointDouble();
        std::cout << "2G = ";
        twoG.print();
        std::cout << "2G on curve: " << (twoG.isOnCurve() ? "YES" : "NO") << std::endl;
        std::cout << std::endl;
    }

    std::cout << "\n========== Modular arithmetic tests ==========\n\n";
    {
        std::cout << "mod(-3, 11) = " << mod(-3LL, 11LL) << " (expected 8)" << std::endl;
        std::cout << "mod(15, 11) = " << mod(15LL, 11LL) << " (expected 4)" << std::endl;
        std::cout << "modInverse(3, 11) = " << modInverse(3LL, 11LL) << " (expected 4)" << std::endl;
        std::cout << "modInverse(7, 11) = " << modInverse(7LL, 11LL) << " (expected 8)" << std::endl;
        std::cout << "modPow(2, 10, 11) = " << modPow(2LL, 10LL, 11LL) << " (expected 1)" << std::endl;
        std::cout << "modPow(3, 5, 11) = " << modPow(3LL, 5LL, 11LL) << " (expected 1)" << std::endl;
    }

    return 0;
}
