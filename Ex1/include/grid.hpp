#pragma once
#include <vector>
#include <string>
#include "sampler.hpp"

enum class AssignmentScheme {
    NGP,  // Nearest Grid Point
    CIC   // Cloud-In-Cell
};

struct RadialBin {
    double r_center;  // Center of the radial bin
    double density;   // Density in the radial bin
    double count;     // Number of particles in the radial bin
};

class Grid3D {
    public:
        // box size: full length along one axis (e.g 500 kpc for [-250, 250])
        // n_cell: 100 or 200
        Grid3D(double box_size, size_t n_cells);

        void reset();

        // task 3: Map particles to grid using the specified assignment scheme
        void assign_mass(const std::vector<Particle>& particles, double particle_mass, AssignmentScheme scheme);

        // Convert cell masses to densities: rho = mass / (dx^3)
        void compute_density();

        // task 4: Spherically average grid densities into radial bins
        std::vector<RadialBin> compute_radial_profile(size_t num_bins, double r_min, double r_max) const;

        // export radial profile to CSV file
        void export_radial_profile(const std::string& filename, const std::vector<RadialBin>& profile) const;

        // helper: cell linear index from (i, j, k)
        inline size_t cell_index(size_t i, size_t j, size_t k) const {
            return i * n_cells_ * n_cells_ + j * n_cells_ + k;
        }

    private:
        double box_size_;  // Full length of the box along one axis
        size_t n_cells_;   // Number of cells along one axis
        double dx_;        // Cell size (box_size / n_cells)
        double cell_volume_; // Volume of a single cell (dx^3)
        double half_box_; // Half the box size (for centering)

        std::vector<double> grid_;  // 3D grid data
};