#pragma once

// Elliptic curve in short Weierstrass form: y^2 = x^3 + ax + b (mod p)
// T = mpz_class for Baby-JubJub, T = long long for simplified cases

// Helper: convert T -> long long (for values known to be small)
inline long long toLongLong(long long v) { return v; }
inline long long toLongLong(const mpz_class& v) { return v.get_si(); }

// Helper: convert long long -> T (avoids ambiguous mpz_class(long long) constructor)
template<typename T> inline T fromLL(long long v);
template<> inline long long fromLL<long long>(long long v) { return v; }
template<> inline mpz_class fromLL<mpz_class>(long long v) {
    if (v >= 0) return mpz_class(static_cast<unsigned long>(v));
    return -mpz_class(static_cast<unsigned long>(-v));
}

template <typename T>
class EllipticCurve {
public:
    T a, b, p;
    T n;           // curve order (#E(F_p), including O_E)
    bool verbose;  // step-by-step output (for long long demonstrations)

    EllipticCurve(const T& a, const T& b, const T& p, const T& n, bool verbose = false)
        : a(a), b(b), p(p), n(n), verbose(verbose) {}

    static std::string toString(const T& val) {
        std::ostringstream oss;
        oss << val;
        return oss.str();
    }

    // Check that the curve is non-singular: Delta = -16(4a^3 + 27b^2) != 0 (mod p).
    // Since -16 is invertible in F_p for any prime p >= 3 (no zero divisors in a field),
    // -16 * X = 0  iff  X = 0, so it suffices to check just 4a^3 + 27b^2 != 0.
    bool isNonSingular() const {
        T disc_core = mod(T(4) * a * a * a + T(27) * b * b, p);
        if (verbose) {
            std::cout << "[isNonSingular] Delta = -16*(4a^3 + 27b^2): 4*" << toString(a) << "^3 + 27*" << toString(b) << "^2 = "
                      << toString(disc_core) << " (mod " << toString(p) << ") "
                      << (disc_core != T(0) ? "!= 0 -> non-singular" : "== 0 -> SINGULAR!")
                      << std::endl;
        }
        return disc_core != T(0);
    }

    // Print the curve equation
    void print() const {
        std::cout << "Elliptic Curve: y^2 = x^3";
        if (a != T(0)) std::cout << " + " << toString(a) << "*x";
        if (b != T(0)) std::cout << " + " << toString(b);
        std::cout << "  (mod " << toString(p) << ")" << std::endl;
    }

    // ===== Bruteforce point enumeration (small p only, T = long long) =====
    // Finds ALL points on the curve by iterating x from 0 to p-1.
    // Stores the order in n. Returns a vector of affine points (excluding O_E).
    // sqrtExampleX: show step-by-step sqrt computation for this x (-1 = don't show).
    std::vector<std::pair<long long, long long>> findAllPointsBruteforce(bool showTable = true, long long sqrtExampleX = 0) {
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
                long long y = modSqrt(rhs[x], pp, false);
                if (y >= 0) {
                    y1arr[x] = y;
                    points.push_back({x, y});
                    points.push_back({x, mod(-y, pp)});
                    count += 2;
                }
            }
        }
        n = static_cast<T>(count);

        // Step-by-step sqrt computation for the chosen x
        if (verbose && sqrtExampleX >= 0 && sqrtExampleX < pp) {
            std::cout << "\nSqrt example for x=" << sqrtExampleX
                      << ": RHS = " << rhs[sqrtExampleX] << " (mod " << pp << ")" << std::endl;
            if (rhs[sqrtExampleX] == 0) std::cout << "  RHS = 0 -> y = 0 (single point)" << std::endl;
            else {
                modSqrt(rhs[sqrtExampleX], pp, true);
                if (y1arr[sqrtExampleX] > 0)
                    std::cout << "  => y1 = " << y1arr[sqrtExampleX] << ", y2 = " << mod(-y1arr[sqrtExampleX], pp) << std::endl;
            }
            std::cout << std::endl;
        }

        // Table of affine coordinates
        if (showTable) {
            const int MAX_COLS = 15, LABEL_W = 7;
            std::cout << "--- All " << count << " affine points of y^2 = x^3 + "
                      << toString(a) << "*x + " << toString(b) << " (mod " << toString(p) << ") ---\n";
            std::cout << "(+ O_E, which has no affine representation)\n\n";
            for (long long start = 0; start < pp; start += MAX_COLS) {
                long long end = std::min(start + MAX_COLS, pp);
                int ncols = static_cast<int>(end - start);
                int sepLen = LABEL_W + 1 + 8 * ncols;
                auto sep = [&]() { std::cout << "|"; for (int i = 0; i < sepLen - 1; i++) std::cout << "-"; std::cout << "|\n"; };
                sep();
                std::cout << "|" << std::setw(LABEL_W) << "x" << "|";
                for (long long x = start; x < end; x++) std::cout << std::setw(6) << x << " |";
                std::cout << "\n"; sep();
                std::cout << "|" << std::setw(LABEL_W) << "RHS" << "|";
                for (long long x = start; x < end; x++) std::cout << std::setw(6) << rhs[x] << " |";
                std::cout << "\n"; sep();
                std::cout << "|" << std::setw(LABEL_W) << "y1" << "|";
                for (long long x = start; x < end; x++) { if (y1arr[x]>=0) std::cout<<std::setw(6)<<y1arr[x]<<" |"; else std::cout<<std::setw(6)<<"-"<<" |"; }
                std::cout << "\n"; sep();
                std::cout << "|" << std::setw(LABEL_W) << "y2" << "|";
                for (long long x = start; x < end; x++) { if (y1arr[x]>0) std::cout<<std::setw(6)<<mod(-y1arr[x],pp)<<" |"; else std::cout<<std::setw(6)<<"-"<<" |"; }
                std::cout << "\n"; sep(); std::cout << "\n";
            }
        }
        return points;
    }

    // ======================== Schoof's Algorithm (templated) ========================
    // Computes #E(F_p) = p + 1 - t, where |t| <= 2*sqrt(p) (Hasse bound).
    // Works for any T: long long (small p) and mpz_class (large p, e.g. Baby-JubJub).
    //
    // Division polynomials hat_psi[n] (stored as dp[n]):
    //   n odd  -> dp[n] = psi_n(x)          (polynomial in x only)
    //   n even -> dp[n] = psi_n(x) / (2y)   (polynomial in x only)

    using PT = PolyT<T>;  // polynomial type with coefficients of type T

    // MSVC C26451: "int+int -> size_t overflow" — false positive.
    // All polynomial indices are bounded by max(deg(psi_l)) = (l^2-1)/2.
    // For Schoof on Baby-JubJub (p ~ 2^254), l < 100, so indices < 5000 — far from INT_MAX.
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 26451)
#endif

    // Compute division polynomials dp[0..maxN] using recurrence relations.
    // dp[0] = 0, dp[1] = 1, dp[2] = 1, dp[3] = 3x^4 + 6ax^2 + 12bx - a^2, ...
    // For k >= 5: recurrence based on k = 2m+1 (odd) or k = 2m (even).
    std::vector<PT> computeDivPolys(int maxN) const {
        T aa = mod(a, p), bb = mod(b, p);
        PT fx = polyNorm(PT{bb, aa, T(0), T(1)});       // f(x) = x^3 + ax + b
        PT fx2 = polyMul(fx, fx, p);                     // f(x)^2
        PT fx2_16 = polyScale(fx2, T(16), p);            // 16 * f(x)^2

        std::vector<PT> dp(maxN + 1, PT{T(0)});
        dp[1] = {T(1)};
        if (maxN >= 2) dp[2] = {T(1)};
        if (maxN >= 3) {
            T a2 = mod(aa * aa, p);
            dp[3] = polyNorm(PT{mod(-a2, p), mod(T(12) * bb, p), mod(T(6) * aa, p), T(0), T(3)});
        }
        if (maxN >= 4) {
            T a2 = mod(aa * aa, p), a3 = mod(a2 * aa, p);
            T b2 = mod(bb * bb, p), ab = mod(aa * bb, p);
            // dp[4] = 2x^6 + 10ax^4 + 40bx^3 - 10a^2*x^2 - 8abx - 2a^3 - 16b^2
            dp[4] = polyNorm(PT{mod(mod(T(-2) * a3, p) - mod(T(16) * b2, p), p),
                                mod(T(-8) * ab, p), mod(T(-10) * a2, p),
                                mod(T(40) * bb, p), mod(T(10) * aa, p), T(0), T(2)});
        }

        for (int k = 5; k <= maxN; k++) {
            if (k % 2 == 1) {
                int m = (k - 1) / 2;  // k = 2m+1
                PT pm3 = polyMul(polyMul(dp[m], dp[m], p), dp[m], p);       // dp[m]^3
                PT pm1_3 = polyMul(polyMul(dp[m+1], dp[m+1], p), dp[m+1], p); // dp[m+1]^3
                if (m % 2 == 0)
                    // m even: dp[2m+1] = 16*f^2 * dp[m+2] * dp[m]^3 - dp[m-1] * dp[m+1]^3
                    dp[k] = polySub(polyMul(fx2_16, polyMul(dp[m+2], pm3, p), p),
                                    polyMul(dp[m-1], pm1_3, p), p);
                else
                    // m odd:  dp[2m+1] = dp[m+2] * dp[m]^3 - 16*f^2 * dp[m-1] * dp[m+1]^3
                    dp[k] = polySub(polyMul(dp[m+2], pm3, p),
                                    polyMul(fx2_16, polyMul(dp[m-1], pm1_3, p), p), p);
            } else {
                int m = k / 2;  // k = 2m
                // dp[2m] = dp[m] * (dp[m+2]*dp[m-1]^2 - dp[m-2]*dp[m+1]^2)
                PT mm1sq = polyMul(dp[m-1], dp[m-1], p);
                PT mp1sq = polyMul(dp[m+1], dp[m+1], p);
                dp[k] = polyMul(dp[m], polySub(polyMul(dp[m+2], mm1sq, p),
                                               polyMul(dp[m-2], mp1sq, p), p), p);
            }
        }
        return dp;
    }

    // X-coordinate of [n]P as num/den (all mod modP, mod p).
    //   n odd:  num = x*dp[n]^2 - 4f*dp[n-1]*dp[n+1],   den = dp[n]^2
    //   n even: num = x*4f*dp[n]^2 - dp[n-1]*dp[n+1],    den = 4f*dp[n]^2
    void xCoordOfScalar(int n, const std::vector<PT>& dp, const PT& fx,
                        const PT& modP, PT& num, PT& den) const {
        PT x_poly = {T(0), T(1)};
        auto R = [&](const PT& a) { return polyMod(a, modP, p); };
        auto M = [&](const PT& a, const PT& b) { return R(polyMul(a, b, p)); };

        if (n == 0) { num = {T(1)}; den = {T(0)}; return; }  // point at infinity

        PT dpn_r = R(dp[n]), dpn1_r = R(dp[n-1]), dpn2_r = R(dp[n+1]);
        PT dpn_sq = M(dpn_r, dpn_r);
        PT adj = M(dpn1_r, dpn2_r);

        if (n % 2 == 1) {
            PT fx4 = polyScale(fx, T(4), p);
            num = R(polySub(M(x_poly, dpn_sq), M(fx4, adj), p));
            den = dpn_sq;
        } else {
            PT fx4 = polyScale(fx, T(4), p);
            PT fx4_dpn_sq = M(fx4, dpn_sq);
            num = R(polySub(M(x_poly, fx4_dpn_sq), adj, p));
            den = fx4_dpn_sq;
        }
    }

    // Y-factor of [n]P: y([n]P) = y * ynum / yden
    //   n = 1:        ynum = 1, yden = 1
    //   n odd  >= 3:  ynum = dp[n+2]*dp[n-1]^2 - dp[n-2]*dp[n+1]^2,  yden = dp[n]^3
    //   n even >= 2:  ynum = dp[n+2]*dp[n-1]^2 - dp[n-2]*dp[n+1]^2,  yden = 16*f^2*dp[n]^3
    void yFactorOfScalar(int n, const std::vector<PT>& dp, const PT& fx,
                         const PT& modP, PT& ynum, PT& yden) const {
        auto R = [&](const PT& a) { return polyMod(a, modP, p); };
        auto M = [&](const PT& a, const PT& b) { return R(polyMul(a, b, p)); };

        if (n <= 1) { ynum = {T(1)}; yden = {T(1)}; return; }

        // Common numerator for both parities
        PT dn1_sq = M(R(dp[n-1]), R(dp[n-1]));
        PT dp1_sq = M(R(dp[n+1]), R(dp[n+1]));
        ynum = R(polySub(M(R(dp[n+2]), dn1_sq), M(R(dp[n-2]), dp1_sq), p));

        // Denominator depends on parity
        PT dpn_r = R(dp[n]);
        PT dpn3 = M(dpn_r, M(dpn_r, dpn_r));   // dp[n]^3

        if (n % 2 == 1) {
            yden = dpn3;
        } else {
            // n even: yden = 16 * f^2 * dp[n]^3
            PT fx_r = R(fx);
            PT fx2 = M(fx_r, fx_r);
            yden = M(polyScale(fx2, T(16), p), dpn3);
        }
    }

    // Main entry point: compute curve order using Schoof's algorithm.
    // Returns #E(F_p) = p + 1 - t.
    T schoofOrder() const {
        PT fx = polyNorm(PT{mod(b, p), mod(a, p), T(0), T(1)});   // f(x) = x^3 + ax + b
        PT x_poly = {T(0), T(1)};                                  // the polynomial x
        T hasse = hasseBound(p);

        // Collect small primes l such that their product > 2 * hasse_bound
        std::vector<long long> primes;
        T prod(1);
        for (long long l = 2; prod <= T(2) * hasse; l++) {
            if (isSmallPrime(l) && fromLL<T>(l) != p) {
                primes.push_back(l);
                prod = prod * fromLL<T>(l);
            }
        }

        std::vector<long long> tvals(primes.size(), 0);

        // ==================== l = 2 ====================
        // t = 0 (mod 2) iff gcd(x^p - x, f(x)) != 1 iff f has a root in F_p
        {
            PT xp = polyPowMod(x_poly, p, fx, p);
            PT g = polyGcd(polySub(xp, x_poly, p), fx, p);
            tvals[0] = (polyDeg(g) > 0) ? 0 : 1;
            if (verbose) std::cout << "[Schoof] l=2: t mod 2 = " << tvals[0] << std::endl;
        }

        // ==================== Odd primes l ====================
        // Schoof equation: pi^2(P) + [q]P = [tau]*pi(P)  in E(F_p[x]/psi_l)
        // where pi(x,y) = (x^p, y^p) is Frobenius, q = p mod l, tau = t mod l.
        //
        // Coordinate representation in the ring F_p[x]/(psi_l):
        //   x-coordinate: rational function num/den (polynomials in x)
        //   y-coordinate: y * (ynum/yden) where ynum, yden are polynomials in x
        //
        // Frobenius y-factors:
        //   y^p   = y * f^{(p-1)/2}      (denoted fph)
        //   y^{p^2} = y * f^{(p^2-1)/2}  (denoted fp2h)

        for (size_t idx = 1; idx < primes.size(); idx++) {
            long long l = primes[idx];
            long long q = toLongLong(mod(p, fromLL<T>(l)));   // q = p mod l (always small)
            int maxDP = (int)(l + 2);
            auto dp = computeDivPolys(maxDP);
            PT psi_l = dp[l];  // l is odd => dp[l] = psi_l

            // Convenience lambdas: R = reduce mod psi_l, M = multiply then reduce
            auto R = [&](const PT& a) { return polyMod(a, psi_l, p); };
            auto M = [&](const PT& a, const PT& b) { return R(polyMul(a, b, p)); };
            auto isZero = [](const PT& a) { return a.size() == 1 && a[0] == T(0); };

            // Frobenius x-coordinates
            PT xp  = polyPowMod(x_poly, p, psi_l, p);     // x^p mod psi_l
            PT xp2 = polyPowMod(xp, p, psi_l, p);         // x^{p^2} mod psi_l

            // Frobenius y-factors
            T fph_exp  = (p - T(1)) / T(2);               // (p-1)/2
            T fp2h_exp = (p * p - T(1)) / T(2);           // (p^2-1)/2
            PT fph  = polyPowMod(fx, fph_exp, psi_l, p);
            PT fp2h = polyPowMod(fx, fp2h_exp, psi_l, p);

            // f(x^p) mod psi_l — needed for x-coordinates when substituting x -> x^p
            PT fp_xp = polySubst(fx, xp, psi_l, p);

            // [q]P coordinates
            PT num_q, den_q, ynum_q, yden_q;
            bool qIsZero = (q == 0);
            if (qIsZero) {
                num_q = {T(0)}; den_q = {T(0)};
                ynum_q = {T(0)}; yden_q = {T(1)};
            } else {
                xCoordOfScalar((int)q, dp, fx, psi_l, num_q, den_q);
                yFactorOfScalar((int)q, dp, fx, psi_l, ynum_q, yden_q);
            }

            bool found = false;
            long long tl = 0;

            // Helper: compute x-coordinate of [tau]*pi(P)
            // Substitutes x^p into division polynomials.
            auto computeTauPiX = [&](int ti, PT& num_t, PT& den_t) {
                auto dpAtXp = [&](int i) -> PT { return polySubst(dp[i], xp, psi_l, p); };
                PT dpti_xp = dpAtXp(ti), dptm_xp = dpAtXp(ti-1), dptp_xp = dpAtXp(ti+1);
                if (ti % 2 == 1) {
                    PT dpti_sq = M(dpti_xp, dpti_xp);
                    PT adj_t = R(polyScale(M(fp_xp, M(dptm_xp, dptp_xp)), T(4), p));
                    num_t = R(polySub(M(xp, dpti_sq), adj_t, p));
                    den_t = dpti_sq;
                } else {
                    PT fx4_xp = polyScale(fp_xp, T(4), p);
                    PT dpti_sq = M(dpti_xp, dpti_xp);
                    PT fx4_sq = M(fx4_xp, dpti_sq);
                    num_t = R(polySub(M(xp, fx4_sq), M(dptm_xp, dptp_xp), p));
                    den_t = fx4_sq;
                }
            };

            // Helper: y-factor of [tau]*pi(P) = f^{(p-1)/2} * ytau_num(x^p) / ytau_den(x^p)
            auto computeTauPiY = [&](int ti, PT& yf_num, PT& yf_den) {
                PT ytn, ytd;
                yFactorOfScalar(ti, dp, fx, psi_l, ytn, ytd);
                ytn = polySubst(ytn, xp, psi_l, p);
                ytd = polySubst(ytd, xp, psi_l, p);
                yf_num = M(fph, ytn);   // f^{(p-1)/2} * ytau_num(x^p)
                yf_den = ytd;            // ytau_den(x^p)
            };

            // ============ Case q = 0: pi^2(P) = [tau]*pi(P) ============
            if (qIsZero) {
                for (long long tau = 1; tau <= (l-1)/2 && !found; tau++) {
                    int ti = (int)tau;
                    PT num_t, den_t;
                    computeTauPiX(ti, num_t, den_t);
                    // Check x: x^{p^2} = num_t/den_t  <=>  xp2*den_t - num_t = 0
                    PT chk = polyNorm(R(polySub(M(xp2, den_t), num_t, p)));
                    if (!isZero(chk)) continue;
                    // x matches! Determine y-sign.
                    PT yf_num, yf_den;
                    computeTauPiY(ti, yf_num, yf_den);
                    // Equal:    fp2h*yf_den - yf_num = 0  =>  tau
                    // Opposite: fp2h*yf_den + yf_num = 0  =>  l - tau
                    PT pos = polyNorm(R(polySub(M(fp2h, yf_den), yf_num, p)));
                    if (isZero(pos)) { tl = tau; found = true; }
                    else {
                        PT neg = polyNorm(R(polyAdd(M(fp2h, yf_den), yf_num, p)));
                        if (isZero(neg)) { tl = l - tau; found = true; }
                    }
                }
                if (!found) tl = 0;  // tau = 0 fallback
            }
            // ============ Case q != 0: pi^2(P) + [q]P = [tau]*pi(P) ============
            else {
                // --- Step 1: check tau = 0 ---
                // tau = 0 => pi^2(P) = -[q]P
                // x: xp2 = num_q/den_q  <=>  xp2*den_q - num_q = 0
                PT x_chk_0 = polyNorm(R(polySub(M(xp2, den_q), num_q, p)));
                if (isZero(x_chk_0)) {
                    // y: pi^2(P) = -[q]P => fp2h = -ynum_q/yden_q
                    //    <=> fp2h*yden_q + ynum_q = 0
                    PT ychk_neg = polyNorm(R(polyAdd(M(fp2h, yden_q), ynum_q, p)));
                    if (isZero(ychk_neg)) { tl = 0; found = true; }
                    else {
                        // pi^2(P) = +[q]P => eigenvalue case: t = +-2*sqrt(q) (mod l)
                        PT ychk_pos = polyNorm(R(polySub(M(fp2h, yden_q), ynum_q, p)));
                        if (isZero(ychk_pos)) {
                            long long w = modSqrt(q % l, l);
                            if (w >= 0 && w != 0) {
                                // Determine sign: is pi(P) = [w]P or [−w]P?
                                PT num_w, den_w;
                                xCoordOfScalar((int)w, dp, fx, psi_l, num_w, den_w);
                                PT wchk = polyNorm(R(polySub(M(xp, den_w), num_w, p)));
                                if (isZero(wchk)) {
                                    PT yw_num, yw_den;
                                    yFactorOfScalar((int)w, dp, fx, psi_l, yw_num, yw_den);
                                    PT yw_chk = polyNorm(R(polySub(M(fph, yw_den), yw_num, p)));
                                    tl = isZero(yw_chk) ? mod(2*w, l) : mod(-2*w, l);
                                } else {
                                    tl = mod(-2*w, l);
                                }
                                found = true;
                            }
                        }
                    }
                }

                // --- Step 2: iterate tau = 1, ..., (l-1)/2 ---
                if (!found) {
                    for (long long tau = 1; tau <= (l-1)/2 && !found; tau++) {
                        int ti = (int)tau;
                        PT num_t, den_t;
                        computeTauPiX(ti, num_t, den_t);

                        // Point addition pi^2(P) + [q]P via chord formula:
                        // DX = x(pi^2 P) - x([q]P),  DY = Y(pi^2 P) - Y([q]P)
                        PT DX = R(polySub(M(xp2, den_q), num_q, p));
                        if (isZero(DX)) continue;  // already handled in step 1
                        PT DY = R(polySub(M(fp2h, yden_q), ynum_q, p));

                        // x3 = f*DY^2*den_q^3 - (xp2*den_q+num_q)*yqd^2*DX^2
                        // (all cross-multiplied to avoid rational function division)
                        PT DY2 = M(DY,DY), DX2 = M(DX,DX);
                        PT yqd2 = M(yden_q,yden_q);
                        PT den_q2 = M(den_q,den_q), den_q3 = M(den_q2,den_q);
                        PT fx_r = R(fx);

                        PT x3n = R(polySub(M(fx_r, M(DY2, den_q3)),
                                           M(M(polyAdd(M(xp2,den_q),num_q,p),yqd2),DX2), p));
                        PT x3d = M(yqd2, M(DX2, den_q));

                        // Compare x: x3n*den_t = num_t*x3d (mod psi_l)?
                        PT cmpL = R(M(x3n, den_t));
                        PT cmpR = R(M(num_t, x3d));
                        PT diff = polyNorm(R(polySub(cmpL, cmpR, p)));
                        if (!isZero(diff)) continue;

                        // x matches! Determine y-sign to distinguish tau vs l-tau.
                        // y3/y = [DY*den_q*(xA - x3) - fp2h*yden_q*DX*x3d] / [yden_q*DX*x3d]
                        PT xA_minus_x3 = R(polySub(M(xp2, x3d), x3n, p));
                        PT y3_num = R(polySub(M(M(DY, den_q), xA_minus_x3),
                                              M(M(fp2h, yden_q), M(DX, x3d)), p));
                        PT y3_den = M(yden_q, M(DX, x3d));

                        // y-factor of [tau]*pi(P)
                        PT yf_num, yf_den;
                        computeTauPiY(ti, yf_num, yf_den);

                        // Compare: y3_num*yf_den vs yf_num*y3_den
                        PT lhs = R(M(y3_num, yf_den));
                        PT rhs_val = R(M(yf_num, y3_den));

                        PT pos_chk = polyNorm(R(polySub(lhs, rhs_val, p)));
                        if (isZero(pos_chk)) { tl = tau; found = true; }
                        else {
                            PT neg_chk = polyNorm(R(polyAdd(lhs, rhs_val, p)));
                            if (isZero(neg_chk)) { tl = l - tau; found = true; }
                        }
                    }
                }
                if (!found) tl = 0;
            }

            tvals[idx] = tl;
            if (verbose) std::cout << "[Schoof] l=" << l << ": t mod " << l << " = " << tl << std::endl;
        }

        // ==================== CRT (Chinese Remainder Theorem) ====================
        T M_crt(1);
        for (auto l : primes) M_crt = M_crt * fromLL<T>(l);

        T t_crt(0);
        for (size_t i = 0; i < primes.size(); i++) {
            T Mi = M_crt / fromLL<T>(primes[i]);
            long long Mi_mod_li = toLongLong(mod(Mi, fromLL<T>(primes[i])));
            long long yi = modInverse(Mi_mod_li, primes[i]);
            t_crt = mod(t_crt + mod(fromLL<T>(tvals[i]) * mod(Mi * fromLL<T>(yi), M_crt), M_crt), M_crt);
        }
        // Shift t_crt from [0, M) to [-M/2, M/2]
        if (t_crt > M_crt / T(2)) t_crt = t_crt - M_crt;

        T order = p + T(1) - t_crt;
        if (verbose) std::cout << "[Schoof] t = " << toString(t_crt) << ", #E = p+1-t = " << toString(order) << std::endl;
        return order;
    }

#ifdef _MSC_VER
#pragma warning(pop)
#endif

};
