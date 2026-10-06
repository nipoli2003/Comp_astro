#pragma once
#include <vector>

struct Particle {
    double x, y, z;
};

class Sampler {
    public:
        // Generates N particles distributed according to rho(r) between r_min and R_max
        static std::vector<Particle> generate_particles(size_t num_particles, unsigned int seed = 42);

        // Total analytical mass between r_min and R_max
        static double get_total_mass();

        // Mass per particle: m_p = M_tot / N
        static double get_particle_mass(size_t num_particles);

    private:
        // Invert the CDF: solve F(r) = u using numerical root finding (bisection or Newton-Raphson)
        static double sample_radius(double u, double total_mass);
};