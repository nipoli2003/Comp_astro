#include "sampler.hpp"
#include "constants.hpp"
#include <random>
#include <cmath>

double Sampler::get_total_mass() {
    return Constants::cumulative_mass(Constants::R_max);
}

double Sampler::get_particle_mass(size_t num_particles) {
    return get_total_mass() / static_cast<double>(num_particles);
}

double Sampler::sample_radius(double u, double total_mass) {
    // Target mass: M(r) = u * M_tot
    double target_m = u * total_mass;

    // Bisection search between r_min and R_max
    double low = Constants::r_min;
    double high = Constants::R_max;
    
    // 35 iterations achieves ~1e-11 precision on [0.1, 250]
    for (int iter = 0; iter < 35; ++iter) {
        double mid = 0.5 * (low + high);
        double m_mid = Constants::cumulative_mass(mid);
        if (m_mid < target_m) low = mid; 
        else high = mid;
    }
    return 0.5 * (low + high);
}

std::vector<Particle> Sampler::generate_particles(size_t num_particles, unsigned int seed) {
    std::vector<Particle> particles(num_particles);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist_u(0.0, 1.0);

    double total_mass = get_total_mass();

    for(size_t i = 0; i < num_particles; i++){
        // 1. Sample radial distance
        double u_r = dist_u(rng);
        double r = sample_radius(u_r, total_mass);

        // 2. Sample isotropic angles
        // cos(theta) is uniform in [-1, 1], phi is uniform in [0, 2*pi)
        double cos_theta = 2.0 * dist_u(rng) - 1.0;
        double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        double phi = 2.0 * M_PI * dist_u(rng);

        // 3. Cartesian positions
        particles[i].x = r * sin_theta * std::cos(phi);
        particles[i].y = r * sin_theta * std::sin(phi);
        particles[i].z = r * cos_theta;
    }

    return particles;
}