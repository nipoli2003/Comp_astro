#include <iostream>
#include <filesystem>
#include "constants.hpp"
#include "sampler.hpp"
#include "grid.hpp"

namespace fs = std::filesystem;

int main() {
    fs::create_directories("output");

    const size_t num_particles = 1'000'000;
    const double box_size = 2.0 * Constants::R_max; // 500 kpc
    const size_t num_radial_bins = 40;

    std::cout << "=== Galaxy Density Profile Sampling ===" << std::endl;
    std::cout << "Central Density rho_0: " << Constants::rho_0 << " M_sun / kpc^3" << std::endl;
    std::cout << "Total mass: " << Sampler::get_total_mass() << " M_sun" << std::endl;
    std::cout << "Generating " << num_particles << " particles..." << std::endl;

    auto particles = Sampler::generate_particles(num_particles);
    double m_p = Sampler::get_particle_mass(num_particles);
    std::cout << "Particle mass: " << m_p << " M_sun" << std::endl;

    // --- SANITY CHECK ON TASK 1 ---
    std::cout << "Running particle sanity check..." << std::endl;
    auto particle_profile = Sampler::compute_particle_profile(
        particles, m_p, num_radial_bins, Constants::r_min, Constants::R_max
    );
    Sampler::export_particle_profile_to_csv("output/profile_particles_raw.csv", particle_profile);
    std::cout << "Exported particle sanity check to output/profile_particles_raw.csv" << std::endl;

    // Grid resolutions to test: 100^3 and 200^3
    std::vector<size_t> grid_sizes = {100, 200};

    for (size_t n_cells : grid_sizes) {
        std::cout << "\nProcessing Grid " << n_cells << "^3 (dx = " << box_size / n_cells << " kpc)..." << std::endl;
        Grid3D grid(box_size, n_cells);

        // 1. Nearest Grid Point (NGP)
        std::cout << "  Running NGP assignment..." << std::endl;
        grid.assign_mass(particles, m_p, AssignmentScheme::NGP);
        grid.compute_density();
        auto profile_ngp = grid.compute_radial_profile(num_radial_bins, Constants::r_min, Constants::R_max);
        
        std::string fn_ngp = "output/profile_ngp_" + std::to_string(n_cells) + ".csv";
        grid.export_profile_to_csv(fn_ngp, profile_ngp);
        std::cout << "  Exported to " << fn_ngp << std::endl;

        // 2. Cloud In Cell (CIC)
        std::cout << "  Running CIC assignment..." << std::endl;
        grid.assign_mass(particles, m_p, AssignmentScheme::CIC);
        grid.compute_density();
        auto profile_cic = grid.compute_radial_profile(num_radial_bins, Constants::r_min, Constants::R_max);
        
        std::string fn_cic = "output/profile_cic_" + std::to_string(n_cells) + ".csv";
        grid.export_profile_to_csv(fn_cic, profile_cic);
        std::cout << "  Exported to " << fn_cic << std::endl;
    }

    std::cout << "\nDone! Run 'python3 plot_results.py' to visualize." << std::endl;
    return 0;
}