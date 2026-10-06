#include "sampler.hpp"
#include "constants.hpp"
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>

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

std::vector<Sampler::ParticleRadialBin> Sampler::compute_particle_profile(
    const std::vector<Particle>& particles,
    double particle_mass,
    size_t num_bins,
    double r_min,
    double r_max) {

    std::vector<ParticleRadialBin> profile(num_bins);
    std::vector<size_t> counts(num_bins, 0);

    double log_min = std::log10(r_min);
    double log_max = std::log10(r_max);
    double d_log_r = (log_max - log_min) / static_cast<double>(num_bins);

    // 1. Bin edges and geometry
    std::vector<double> r_edges(num_bins + 1);
    for (size_t b = 0; b <= num_bins; ++b) {
        r_edges[b] = std::pow(10.0, log_min + b * d_log_r);
    }

    for (size_t b = 0; b < num_bins; ++b) {
        profile[b].r_low = r_edges[b];
        profile[b].r_high = r_edges[b + 1];
        profile[b].r_center = std::sqrt(r_edges[b] * r_edges[b + 1]);

        // Exact mass enclosed in this shell
        double dM_exact = Constants::cumulative_mass(r_edges[b + 1]) - Constants::cumulative_mass(r_edges[b]);
        double shell_volume = (4.0 / 3.0) * M_PI * (std::pow(r_edges[b + 1], 3) - std::pow(r_edges[b], 3));
        
        // Exact average density across the shell volume
        profile[b].rho_exact_shell = dM_exact / shell_volume;
    }

    // 2. Count particles in shells
    for (const auto& p : particles) {
        double r = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        if (r >= r_min && r < r_max) {
            size_t b = static_cast<size_t>((std::log10(r) - log_min) / d_log_r);
            if (b < num_bins) {
                counts[b]++;
            }
        }
    }

    // 3. Assign measured densities
    for (size_t b = 0; b < num_bins; ++b) {
        double r1 = profile[b].r_low;
        double r2 = profile[b].r_high;
        double shell_volume = (4.0 / 3.0) * M_PI * (r2 * r2 * r2 - r1 * r1 * r1);

        profile[b].count = counts[b];
        profile[b].density = (counts[b] * particle_mass) / shell_volume;
    }

    return profile;
}

void Sampler::export_particle_profile_to_csv(
    const std::string& filename,
    const std::vector<ParticleRadialBin>& profile) {

    std::ofstream out(filename);
    out << std::scientific << std::setprecision(8);
    out << "r_low,r_high,r_center,rho_measured,rho_exact_shell,count\n";
    for (const auto& bin : profile) {
        if (bin.count > 0) {
            out << bin.r_low << ","
                << bin.r_high << ","
                << bin.r_center << ","
                << bin.density << ","
                << bin.rho_exact_shell << ","
                << bin.count << "\n";
        }
    }
}