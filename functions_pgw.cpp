#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include "bessel_funcs.h"
// #include "trigonometric_funcs.h"

constexpr int N = 101; // Must be odd for Simpson
constexpr double PI = M_PI;
const double eps = 0.001;


double Greens_functions(double eta1, double eta2, double k, double cs = 1.0 / std::sqrt(3.0)) {
    // This function evaluates G' * G' in the limit k*eta >> 1
    // Parameters: - eta1: first integrated time variable
    //             - eta2: second integrated time variable
    //             - k: gravitational wave wavenumber
    //             - cs: speed of sound
    double Greens;
    double nu = (1.0 - 3.0 * cs * cs) / (1.0 + 3.0 * cs * cs);
    if (cs >= 1.0 / std::sqrt(3.0) - eps) {
        Greens = 0.5 * std::cos(k * (eta1 - eta2));
    } else {
        Greens = 0.5 * (k*eta1) * (k*eta2) * (spherical_bessel_y(nu, k*eta1)*spherical_bessel_y(nu, k*eta2) + spherical_bessel_j_adaptive(nu, k*eta1)*spherical_bessel_j_adaptive(nu, k*eta2) );
    }

    return Greens;
}



double Pv_tilde(double x) {
    // Analytic approximation for the spectral density of sound wave plane wave amplitudes
    // Parameter: - x: sound wave dimensionless wavenumber (x = p R_*)
    double xr = x / (2.0 * PI);
    return (3.0 * PI) / std::pow(2.0 * PI, 3) * xr * xr / (1.0 + std::pow(xr, 6));
}



double kernel_integrand_mu(double tau1, double tau2, double mu, double x, double z,
                           double tau_star, double cs = 1.0 / std::sqrt(3.0)) {
    // Integrand of the kernel function (see equation 2.33 of arXiv:2409.01426 [gr-qc])
    // Parameters: - tau1: first integrated dimensionless time variable (tau1 = eta1/R_*)
    //             - tau2: second integrated dimensionless time variable (tau2 = eta2/R_*)
    //             - mu: direction cosine (mu = \hat{p} \cdot \hat{k}) 
    //             - x: sound wave dimensionless wavenumber (x = p R_*)
    //             - z: dimensionless gravitational wave wavenumber (z = k R_*)
    //             - tau_star: dimensionless starting time of the acoustic source (tau_star = eta_star / R_*)
    //             - cs: speed of sound
    double y = std::sqrt(x * x + z * z - 2.0 * mu * x * z);
    double nu = (1.0 - 3.0 * cs * cs) / (1.0 + 3.0 * cs * cs);
    double tau_minus = tau1 - tau2;
    return std::pow(tau_star, 2 * nu) * std::pow(tau1 * tau2, -1.0 - nu) *
           Greens_functions(tau1, tau2, z, cs) *
           std::cos(cs * x * tau_minus) * std::cos(cs * y * tau_minus);
}



double integrand_Pgw_spec_dens_mu_interpolation(double mu, double x, double tau1, double tau2,
                                   double z, double tau_star, double tau_end, double Pv_x, double Pv_y,
                                   double cs = 1.0 / std::sqrt(3.0)) {
    // Integrand of the spectral density of tensor modes (see equation 2.36 of arXiv:2409.01426 [gr-qc])
    // Parameters: - tau1: first integrated dimensionless time variable (tau1 = eta1/R_*)
    //             - tau2: second integrated dimensionless time variable (tau2 = eta2/R_*)
    //             - mu: direction cosine (mu = \hat{p} \cdot \hat{k}) 
    //             - x: sound wave dimensionless wavenumber (x = p R_*)
    //             - z: dimensionless gravitational wave wavenumber (z = k R_*)
    //             - tau_star: dimensionless starting time of the acoustic source (tau_star = eta_star / R_*)
    //             - tau_end: dimensionless final time of the acoustic source (tau_end = eta_star / R_*)
    //             - Pv_x: spectral density of sound wave amplitudes evaluated in x
    //             - Pv_y: spectral density of sound wave amplitudes evaluated in y = sqrt(x^2 +z^2 - 2mu*x*z)
    //             - cs: speed of sound
    // 
    // The end of the acoustic phase here is arbitrary, but it should be consistent with the lifetime of the source $H_* tau_v = (H_* R_*)/\bar{U} 
    double nu = (1.0 - 3.0 * cs * cs) / (1.0 + 3.0 * cs * cs);
    double y = std::sqrt(x * x + z * z - 2.0 * mu * x * z);
    double rho = std::pow(1.0 - mu * mu, 2) * std::pow(x, 4) * std::pow(z, 3)  / (y * y);

    // double y_r = y / (2.0 * PI);
    // double x_r = x / (2.0 * PI);
    // double Pv_x = 3.0 * PI * std::pow(2.0 * PI, -3) * x_r * x_r / (1.0 + std::pow(x_r, 6));
    // double Pv_y = 3.0 * PI * std::pow(2.0 * PI, -3) * y_r * y_r / (1.0 + std::pow(y_r, 6));

    return tau_star / (std::pow(PI, 2) * std::pow(z, 3)) * Pv_x * Pv_y * rho *
           kernel_integrand_mu(tau1, tau2, mu, x, z, tau_star, cs);
}



// Define linear interpolation given x_data and y_data
double interpolate(const std::vector<double>& x_data, const std::vector<double>& y_data, double x) {
    // Obtains the value of the function y(x) from the linear interpolation between a table of given value x_data and y_data
    // Parameters: - x_data: grid of x values for the interpolation
    //             - y_data: grid of y values for the interpolation
    //             - x: position at which we need to evaluate the function y(x)
    // Returns: y(x)
    // 
    // I use this function to evalute the spectral density of sound wave amplitudes Pv(x). 
    // I use PTtools to create a .csv file that contains the values of Pv(x) for x in [0.001, 1000].
    // In the main function I read this .csv file, which forms the grid x_data and y_data. 
    // In the integration "power_spectrum_mu_4D_interpolation" I finally use this "interpolate" function to get the value of Pv(x) for arbitrary x
    if (x <= x_data.front()) return y_data.front();
    if (x >= x_data.back())  return y_data.back();

    auto it = std::upper_bound(x_data.begin(), x_data.end(), x);
    size_t i = std::distance(x_data.begin(), it) - 1;

    double x0 = x_data[i], x1 = x_data[i+1];
    double y0 = y_data[i], y1 = y_data[i+1];

    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

// === Simpson's Rule ===

double simpson_1d(const std::vector<double>& y, const std::vector<double>& x) {
    size_t N = y.size();
    if (N < 3 || N % 2 == 0)
        throw std::invalid_argument("Simpson's rule requires an odd number of points >= 3.");

    double result = 0.0;
    for (size_t i = 0; i < N - 2; i += 2) {
        double h = x[i + 2] - x[i];
        result += h / 6.0 * (y[i] + 4.0 * y[i + 1] + y[i + 2]);
    }
    return result;
}





std::vector<double> power_spectrum_mu_4D_interpolation(double kp, double eta_end, double omega, 
                                                        const std::vector<double>& k, const std::vector<double>& mu_vals, 
                                                        const std::vector<double>& x_vals, const std::vector<double>& tau1_vals, 
                                                        const std::vector<double>& tau2_vals, const std::vector<double>& x_data, 
                                                        const std::vector<double>& Pv_data) {
    // Spectral density of tensor modes (see equation 2.36 of arXiv:2409.01426 [gr-qc])
    // Parameters: - kp: peak amplitude (kp = 2 pi/ R_*)
    //             - eta_end: final time of the acoustic source (eta_star = 1)
    //             - omega: EOS parameter (p = omega e)
    //             - k: vector of gravitational wave wavenumbers to sample
    //             - mu_vals: vector of mu-values that form the grid for simpson integration (the larger the number of element to more precise the integration)
    //             - x_vals: vector of x-values that form the grid for simpson integration (the larger the number of element to more precise the integration)
    //             - tau1_vals: vector of tau1-values that form the grid for simpson integration (the larger the number of element to more precise the integration)
    //             - tau2_vals: vector of tau2-values that form the grid for simpson integration (the larger the number of element to more precise the integration)
    //             - x_data: vector of sound wave dimensionless wavenumbers at which I sample the spectral density Pv (obtained from PTtools)
    //             - Pv_data: vector of values of the spectral density of sound wave amplitudes Pv corresponding to x_data
    // The end of the acoustic phase (eta_end) here is arbitrary, but it should be consistent with the lifetime of the source $H_* tau_v = (H_* R_*)/\bar{U} 
   
    // Initialize the grid for simpson integration
    size_t N_mu = mu_vals.size();
    size_t N_x = x_vals.size();
    size_t N_t1 = tau1_vals.size();
    size_t N_t2 = tau2_vals.size();
    // Initialize an empty vector to take the values of the gravitational wave spectral density
    std::vector<double> PS(k.size(), 0.0);

    double Lf = 2 * M_PI / kp; // length scale of the fluid ~ R_*
    double cs = std::sqrt(omega);
    double tau_star = 1.0 / Lf;            // eta_start = 1
    double tau_end = eta_end / Lf;

    // Start a for-loop over each value of k. For every k, I evaluate the integrand of equation 2.36 (arXiv:2409.01426 [gr-qc]) and perform the 4D integration
    for (size_t i = 0; i < k.size(); ++i) { 
        std::cout << "Processing k[" << i << "] = " << k[i] << std::endl;
        double z = k[i] * Lf;

        // I use 4 nested for-loops to initialize the integrand on the entire grid. 
        // In every loop I set the value of one integration variable (mu, x, tau1, tau2)
        // First loop: mu
        std::vector<double> integrand_mu(N_mu); // 1D vector that contains the value of the integrand for every value of mu on the grid (x, tau1, tau2 integrated out)
        for (size_t j = 0; j < N_mu; ++j) {
            double mu = mu_vals[j]; // fix mu

            // Second loop: x
            std::vector<double> integrand_x(N_x); // 1D vector that contains the value of the integrand for fixed (mu) and every value of x on the grid (tau1, tau2 integrated out)
            for (size_t m = 0; m < N_x; ++m) {
                double x = x_vals[m]; // fix x (and thus y)
                double y = std::sqrt(x * x + z * z - 2.0 * mu * x * z);
                double Pv_x = interpolate(x_data, Pv_data, x); // initialize the value of the spectral density of plane wave amplitudes at x by 
                                                               // interpolating the tabulated values x_data Pv_data obtained with PTtools. 
                                                               // Upload x_data and Pv_data from .csv file in main() 
                double Pv_y = interpolate(x_data, Pv_data, y);

                // Third loop: tau1
                std::vector<double> integrand_tau1(N_t1); // 1D vector that contains the value of the integrand for fixed (mu, x) and every value of tau1 on the grid (tau2 integrated out)
                for (size_t t1 = 0; t1 < N_t1; ++t1) {
                    double tau1 = tau1_vals[t1]; // fix tau1

                    // Fourth loop (innermost): tau2
                    std::vector<double> integrand_tau2(N_t2); // 1D vector that contains the value of the integrand for fixed (mu, x, tau1) and every value of tau2 on the grid
                    for (size_t t2 = 0; t2 < N_t2; ++t2) {
                        double tau2 = tau2_vals[t2]; // fix tau2

                        // Having fixed mu, x, tau1, we now compute the integrand for every value of tau2 on the grid.
                        integrand_tau2[t2] = integrand_Pgw_spec_dens_mu_interpolation(
                            mu, x, tau1, tau2, z,
                            tau_star, tau_end, Pv_x, Pv_y, cs
                        );
                    }
                    // Perform the first integration over tau2
                    integrand_tau1[t1] = simpson_1d(integrand_tau2, tau2_vals); // this is 1D vector of the integrand where (mu, x) are fixed and tau1 has values on the grid
                }
                // Perform the second integration over tau1
                integrand_x[m] = simpson_1d(integrand_tau1, tau1_vals); // this is 1D vector of the integrand where (mu) is fixed and x has values on the grid
            }
            // Perform the third integration over x
            integrand_mu[j] = simpson_1d(integrand_x, x_vals); // this is 1D vector of the integrand where mu has values on the grid
        }
        // Perform the fourth integration over mu
        PS[i] = simpson_1d(integrand_mu, mu_vals);
    }

    return PS; // returns a vector containing the result of the 4D integration for every value of k
}





