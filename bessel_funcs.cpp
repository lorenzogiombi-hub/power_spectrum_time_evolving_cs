// In this document I define the important Bessel functions. 
// Bessel functions can be more easily obtained from <boost/math/special_functions/bessel.hpp>
// If this package is not available (as it is for my laptop), the following definitions can be used

#include <cmath>
#include <vector>
#include <iostream>
// #include <boost/math/special_functions/bessel.hpp>


constexpr double PI = M_PI;

// Define logspace and linspace for convenience
std::vector<double> linspace(double start, double end, int num) {
    // Analogous to numpy.linspace(start, end, numer_of_points)
    std::vector<double> result;
    if (num <= 1) {
        result.push_back(start);
        return result;
    }

    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + step * i);
    }
    return result;
}

std::vector<double> logspace(double start_exp, double end_exp, int num) {
    // Analogous to numpy.logspace(start, end, numer_of_points)
    std::vector<double> result;
    if (num <= 1) {
        result.push_back(std::pow(10.0, start_exp));
        return result;
    }

    double step = (end_exp - start_exp) / (num - 1);
    for (int i = 0; i < num; ++i) {
        double exponent = start_exp + i * step;
        result.push_back(std::pow(10.0, exponent));
    }
    return result;
}


// // Compute Pochhammer symbol (a)_k
// double pochhammer(double a, int k) {
//     double result = 1.0;
//     for (int i = 0; i < k; ++i)
//         result *= a + i;
//     return result;
// }


// // Gamma function via Lanczos or use tgamma from <cmath>
// double gamma(double z) {
//     return std::tgamma(z);  // built-in gamma function
// }


// _0F1(a; z) hypergeometric function
double hypergeometric_0F1(double a, double z, int terms = 50) {
    double sum = 1.0;
    double term = 1.0;
    for (int k = 1; k < terms; ++k) {
        term *= z / (k * (a + k - 1));
        sum += term;
        if (std::fabs(term) < 1e-12) break;  // early stopping
    }
    return sum;
}


// Bessel function J_nu(x) via hypergeometric series ---- from Wikipedia https://en.wikipedia.org/wiki/Bessel_function
// This relation is exact, but hypergeometric_0F1 is numerically stable only for small enough x, 
double bessel_J_hyper(double nu, double x) {
    double coeff = std::pow(x / 2.0, nu) / std::tgamma(nu + 1.0);
    double z = -0.25 * x * x;
    return coeff * hypergeometric_0F1(nu + 1.0, z);
}


// --- Asymptotic Approximation for J_nu (large x) ---  from Wikipedia https://en.wikipedia.org/wiki/Bessel_function
double bessel_J_asymptotic(double nu, double x) {
    double phi = x - (nu * PI / 2.0) - (PI / 4.0);
    return std::sqrt(2.0 / (PI * x)) * std::cos(phi);
}


// --- Adaptive Bessel J_nu ---  
double bessel_J_adaptive(double nu, double x) {
    // Bessel function of the first kind J_nu(x)
    // Use bessel_J_hyper for small enough x and bessel_J_asymptotic for large x
    if (x < 30.0)
        return bessel_J_hyper(nu, x);
    else
        return bessel_J_asymptotic(nu, x);
}



// Bessel Y_nu(x) for non-integer nu using the identity
double bessel_Y_from_J(double nu, double x) {
    // Compute the Bessel function of the second kind Y_nu(x) from the bessel function of the first kind J_nu(x)
    double Jnu = bessel_J_adaptive(nu, x);
    double Jneg_nu = bessel_J_adaptive(-nu, x);
    double sin_nupi = std::sin(nu * M_PI);

    if (std::abs(sin_nupi) < 1e-12) {
        std::cerr << "Y_nu undefined for integer nu: singularity at nu = " << nu << "\n";
        return NAN;
    }

    return (Jnu * std::cos(nu * M_PI) - Jneg_nu) / sin_nupi;
}



// --- Final: Spherical Bessel j_nu ---

double spherical_bessel_j(double nu, double x) {
    // Spherical Bessel function of first kind j_nu(x) for small x
    if (x == 0.0) return (nu == 0.0) ? 1.0 : 0.0;
    return std::sqrt(M_PI / (2.0 * x)) * bessel_J_hyper(nu + 0.5, x); // j_n(x) = sqrt(pi/2x) J_[n+1/2](x)
}

double spherical_bessel_j_approx(double nu, double x) {
    // Spherical Bessel function of first kind j_nu(x) for large x
    if (x == 0.0) return (nu == 0.0) ? 1.0 : 0.0;
    return std::sqrt(M_PI / (2.0 * x)) * bessel_J_asymptotic(nu + 0.5, x); // j_n(x) = sqrt(pi/2x) J_[n+1/2](x)
}




double spherical_bessel_j_adaptive(double nu, double x) {
    // Spherical Bessel function of first kind j_nu(x)
    // Uses the exact relation bessel_J_hyper for small x and the large argument approximation for large x
    if (x == 0.0) return (nu == 0.0) ? 1.0 : 0.0;
    return std::sqrt(PI / (2.0 * x)) * bessel_J_adaptive(nu + 0.5, x); // j_n(x) = sqrt(pi/2x) J_[n+1/2](x)
}


// Spherical Bessel function of the second kind
double spherical_bessel_y(double nu, double x) {
    // Spherical Bessel function of second kind y_nu(x) from J_nu(x)
    if (x == 0.0) return -INFINITY;  // Singularity at x = 0
    return std::sqrt(M_PI / (2.0 * x)) * bessel_Y_from_J(nu + 0.5, x);   // y_n(x) = sqrt(pi/2x) Y_[n+1/2](x)
}
