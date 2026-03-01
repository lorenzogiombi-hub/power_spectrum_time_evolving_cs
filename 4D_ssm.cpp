#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>  // Already included, but make sure this is present
#include <sstream>
// #include <boost/math/special_functions/bessel.hpp>
#include "bessel_funcs.h"
// #include "integrand_funcs.h"
#include "functions_pgw.h"



constexpr int N = 201; // Must be odd for Simpson
constexpr double PI = M_PI;


// Read CSV file and return x and y vectors
bool read_csv(const std::string& filename, std::vector<double>& x_vals, std::vector<double>& y_vals) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file '" << filename << "'\n";
        return false;
    }

    std::string line;
    bool header_skipped = false;
    while (std::getline(file, line)) {
        if (!header_skipped) { header_skipped = true; continue; } // Skip header

        std::stringstream ss(line);
        std::string x_str, y_str;
        if (std::getline(ss, x_str, ',') && std::getline(ss, y_str)) {
            x_vals.push_back(std::stod(x_str));
            y_vals.push_back(std::stod(y_str));
        }
    }

    file.close();
    return true;
}




int main() {

    double omega = 1/3.0;
    double HR = 0.2;
    double k_peak = 2.0 * PI /HR  * 2.0 / (1.0 + 3.0 * omega); // k_peak = 2pi / Lf, where Lf = 2pi/k_peak
    double Lf = 2.0 * PI / k_peak;

    
    // double omega = 0.15;
    double eta_start = 1.0;
    // double Ubar = 0.025304304535680625; // vw = 0.44, alpha = 0.016
    // double Ubar = 0.1379016841186718; // vw = 0.82, alpha = 0.11
    double Ubar = 0.41702219624041087; // vw = 0.56, alpha = 0.40

    double Ht = HR/Ubar; 
    double eta_end = eta_start + Ht;
    // double eta_sh = Lf / (4.0 * PI * std::sqrt(3.0));
    // double Nsh_val = 5000.0;
    // double eta_end = eta_start + Nsh_val * eta_sh;
    // double eta_end = 1.708853815426791; // vw = 0.8, alpha = 0.01, HR = 0.01, Ubar = 0.014107281053399057
    // double eta_end = 1.957964517; // vw = 0.3, alpha = 0.01, HR = 0.01, Ubar = 0.01043881728381361
    // double eta_end = 8.088538154267908;
    double tau_star = eta_start/Lf;
    double tau_end = eta_end/Lf;

    std::vector<double> x_data, Pv_data;
    std::string filename_input = "Velocity_data_cutting/vw_0.56_alpha_0.400.csv"; // Change path as needed
    if (!read_csv(filename_input, x_data, Pv_data)) {
        return 1;
    }

    std::vector<double> k_new = logspace(-1, 3, 80);       // You define logspace
    // std::vector<double> k_new = logspace(1, 4, 50);       // You define logspace
    std::vector<double> mu_vals = linspace(-1, 1, N);
    std::vector<double> x_vals = logspace(-3, 3, N);
    std::vector<double> tau1_vals = linspace(tau_star, tau_end, N);
    std::vector<double> tau2_vals = linspace(tau_star, tau_end, N);

    std::vector<double> Ps = power_spectrum_mu_4D_interpolation(k_peak, eta_end, omega, k_new, mu_vals, x_vals, tau1_vals, tau2_vals, x_data, Pv_data);


    // auto Ps = power_spectrum_mu_4D_new(k_peak, eta_end, omega, k);


    std::cout << "\nReduced_P_gw_list_mu_4D:\n";
    for (size_t i = 0; i < Ps.size(); ++i)
        std::cout << "k[" << i << "] => " << std::setprecision(10) << Ps[i]  << std::endl;
        // std::cout << "k[" << i << "] => " << std::setprecision(10) << P[i] <<  " or " << Ps[i] << std::endl;


    // Save k and PS into a CSV file
    // Format filename
    std::ostringstream filename;
    // filename << "DataSSM/Nsh_1500/data_kp200pi/Pk_N50_omega_" << std::fixed << std::setprecision(2) << omega << ".csv";
    filename << "Pk_HR_0.2_vw_0.56_alpha_0.40_omega_0.333.csv";



    // std::ofstream outfile("power_spectrum_output.csv");
    std::ofstream outfile(filename.str());


    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open output CSV file!" << std::endl;
        return 1;
    }

    // Write header
    outfile << "k,P(k)\n";

    // Write data
    for (size_t i = 0; i < k_new.size(); ++i) {
        outfile << k_new[i] << "," << Ps[i] << "\n";
    }

    outfile.close();
    std::cout << "Results saved to '" << filename.str() << "'\n";

    return 0;
}
