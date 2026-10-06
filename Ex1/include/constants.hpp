#pragma once
#include <cmath>

namespace Constants {
    // Input parameters from sheet
    constexpr double r_0 = 1.0;             // Core radius in kpc
    constexpr double sigma_v = 150.0;       // Velocity scale in km/s
    constexpr double R_max = 250.0;         // Max galaxy radius in kpc
    constexpr double r_min = 0.1 * r_0;     // Min galaxy radius (0.1 kpc)

    // Gravitational constant G in (kpc * (km/s)^2) / M_sun
    constexpr double G = 4.30091727e-6;

    // Derived rho_0: rho_0 = 9 * sigma_v^2 / (4 * pi * G * r_0^2)
    // Yields ~3.746e9 M_sun / kpc^3
    const double rho_0 = (9.0 * sigma_v * sigma_v) / 
                         (4.0 * M_PI * G * r_0 * r_0);

    // Analytical profile rho(r) = rho_0 / (1 + (r/r_0)^2)
    inline double analytical_density(double r) {
        double ratio = r / r_0;
        return rho_0 / (1.0 + ratio * ratio);
    }

    // Cumulative mass function M(r) = Integral(4 * pi * r'^2 * rho(r') dr') from r_min to r
    // Indefinite integral: 4 * pi * rho_0 * r_0^3 * [ (r/r_0) - atan(r/r_0) ]
    inline double indefinite_mass_integral(double r) {
        double x = r / r_0;
        return 4.0 * M_PI * rho_0 * std::pow(r_0, 3) * (x - std::atan(x));
    }

    inline double cumulative_mass(double r) {
        return indefinite_mass_integral(r) - indefinite_mass_integral(r_min);
    }
}