#include "grid.hpp"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

Grid3D::Grid3D(double box_size, size_t n_cells)
    : box_size_(box_size), n_cells_(n_cells) {
    dx_ = box_size_ / static_cast<double>(n_cells_);
    cell_volume_ = dx_ * dx_ * dx_;
    half_box_ = 0.5 * box_size_;
    grid_.resize(n_cells_ * n_cells_ * n_cells_, 0.0);
}

void Grid3D::reset() {
    std::fill(grid_.begin(), grid_.end(), 0.0);
}

void Grid3D::assign_mass(const std::vector<Particle>& particles, double particle_mass, AssignmentScheme scheme) {
    reset();

    for(const auto& p : particles) {
        double shifted_x = p.x + half_box_;
        double shifted_y = p.y + half_box_;
        double shifted_z = p.z + half_box_;

        if(scheme == AssignmentScheme::NGP){
            // Nearest Grid Point: round to cell index
            int i = static_cast<int>(std::floor(shifted_x / dx_));
            int j = static_cast<int>(std::floor(shifted_y / dx_));
            int k = static_cast<int>(std::floor(shifted_z / dx_));

            // Guard boundary checks
            if(i >= 0 && i < static_cast<int>(n_cells_) &&
               j >= 0 && j < static_cast<int>(n_cells_) &&
               k >= 0 && k < static_cast<int>(n_cells_)) {
                grid_[index(i, j, k)] += particle_mass;
            }
        }

        else if(scheme == AssignmentScheme::CIC){
            // Cloud-In-Cell: coordinates relative to cell centers
            double u = shifted_x / dx_ - 0.5;
            double v = shifted_y / dx_ - 0.5;
            double w = shifted_z / dx_ - 0.5;

            int i = static_cast<int>(std::floor(u));
            int j = static_cast<int>(std::floor(v));
            int k = static_cast<int>(std::floor(w));

            // Fractional offsets
            double dx_part = u - i;
            double dy_part = v - j;
            double dz_part = w - k;

            // 1D shape factors
            double tx[2] = {1.0 - dx_part, dx_part};
            double ty[2] = {1.0 - dy_part, dy_part};
            double tz[2] = {1.0 - dz_part, dz_part};

            // Distribute across the 8 neighboring vertices
            for (int di = 0; di < 2; ++di) {
                int gi = i + di;
                if (gi < 0 || gi >= static_cast<int>(n_cells_)) continue;
                for (int dj = 0; dj < 2; ++dj) {
                    int gj = j + dj;
                    if (gj < 0 || gj >= static_cast<int>(n_cells_)) continue;
                    for (int dk = 0; dk < 2; ++dk) {
                        int gk = k + dk;
                        if (gk < 0 || gk >= static_cast<int>(n_cells_)) continue;

                        double weight = tx[di] * ty[dj] * tz[dk];
                        grid_[index(gi, gj, gk)] += particle_mass * weight;
                    }
                }
            }
        }

        else {
            throw std::invalid_argument("Unknown assignment scheme");
        }
    }
}

void Grid3D::compute_density() {
    for(auto& val : grid_) {
        val /= cell_volume_;
    }
}

std::vector<RadialBin> Grid3D::compute_radial_profile(size_t num_bins, double r_min, double r_max) const {
    std::vector<RadialBin> profile(num_bins);
    std::vector<double> sum_density(num_bins, 0.0);
    std::vector<size_t> counts(num_bins, 0);

    double log_min = std::log10(r_min);
    double log_max = std::log10(r_max);
    double d_log_r = (log_max - log_min) / static_cast<double>(num_bins);

    // Initialize bin centers
    for (size_t b = 0; b < num_bins; ++b) {
        double r_low = std::pow(10.0, log_min + b * d_log_r);
        double r_high = std::pow(10.0, log_min + (b + 1) * d_log_r);
        profile[b].r_center = 0.5 * (r_low + r_high);
    }

    // Accumulate each cell center into its spherical radial bin
    for (size_t i = 0; i < n_cells_; ++i) {
        double x = (static_cast<double>(i) + 0.5) * dx_ - half_box_;
        for (size_t j = 0; j < n_cells_; ++j) {
            double y = (static_cast<double>(j) + 0.5) * dx_ - half_box_;
            for (size_t k = 0; k < n_cells_; ++k) {
                double z = (static_cast<double>(k) + 0.5) * dx_ - half_box_;

                double r = std::sqrt(x * x + y * y + z * z);
                if (r >= r_min && r < r_max) {
                    size_t b = static_cast<size_t>((std::log10(r) - log_min) / d_log_r);
                    if (b < num_bins) {
                        sum_density[b] += grid_[index(i, j, k)];
                        counts[b] += 1;
                    }
                }
            }
        }
    }

    for (size_t b = 0; b < num_bins; ++b) {
        profile[b].count = counts[b];
        profile[b].density = (counts[b] > 0) ? (sum_density[b] / counts[b]) : 0.0;
    }

    return profile;
}

void Grid3D::export_profile_to_csv(const std::string& filename, const std::vector<RadialBin>& profile) const {
    std::ofstream out(filename);
    out << std::scientific << std::setprecision(8);
    out << "r,rho,count\n";
    for (const auto& bin : profile) {
        if (bin.count > 0) {
            out << bin.r_center << "," << bin.density << "," << bin.count << "\n";
        }
    }
}