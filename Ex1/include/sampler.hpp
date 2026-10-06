#pragma once
#include <vector>
#include <string>

struct Particle {
    double x, y, z;
};

class Sampler {
    public:
        struct ParticleRadialBin {
            double r_low;
            double r_high;
            double r_center;
            double density;
            double rho_exact_shell; // Exact integrated mass in shell / shell volume
            size_t count;
        };

        // Generates N particles distributed according to rho(r) between r_min and R_max
        static std::vector<Particle> generate_particles(size_t num_particles, unsigned int seed = 42);
        // Total analytical mass between r_min and R_max
        static double get_total_mass();
        // Mass per particle: m_p = M_tot / N
        static double get_particle_mass(size_t num_particles);

        static std::vector<ParticleRadialBin> compute_particle_profile(
            const std::vector<Particle>& particles,
            double particle_mass,
            size_t num_bins,
            double r_min,
            double r_max
        );

        static void export_particle_profile_to_csv(
            const std::string& filename,
            const std::vector<ParticleRadialBin>& profile
        );

    private:
        // Invert the CDF: solve F(r) = u using numerical root finding (bisection or Newton-Raphson)
        static double sample_radius(double u, double total_mass);
};