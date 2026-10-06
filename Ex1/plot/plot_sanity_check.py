import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("output/profile_particles_raw.csv", delimiter=",", names=True)

r = data['r_center']
rho = data['rho_measured']
rho_exact = data['rho_exact_shell']
counts = data['count']

# Poisson relative uncertainty: sigma / rho = 1 / sqrt(N)
rel_error = (rho - rho_exact) / rho_exact
poisson_1sigma = 1.0 / np.sqrt(counts)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True, gridspec_kw={'height_ratios': [2.5, 1]})

# Density Profile
ax1.plot(r, rho_exact, 'k--', label=r'Exact Shell Average $\bar{\rho}_{\mathrm{exact}}$', lw=1.5)
ax1.errorbar(r, rho, yerr=rho * poisson_1sigma, fmt='o-', color='navy', markersize=3, 
             capsize=2, elinewidth=0.8, label=r'Sampled Particles ($N=10^6$)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$\rho(r)\; [M_\odot / \mathrm{kpc}^3]$')
ax1.set_title('Task 1 Sanity Check: Exact Shell Density vs Poisson Errors')
ax1.grid(True, which="both", ls=":", alpha=0.5)
ax1.legend()

# Residuals with shaded 1-sigma and 2-sigma Poisson confidence bands
ax2.plot(r, rel_error, 'o-', color='navy', markersize=3, label=r'Residual $(\rho - \bar{\rho})/\bar{\rho}$')
ax2.fill_between(r, -poisson_1sigma, poisson_1sigma, color='gray', alpha=0.3, label=r'$\pm 1\sigma$ Poisson noise ($1/\sqrt{N}$)')
ax2.fill_between(r, -2 * poisson_1sigma, 2 * poisson_1sigma, color='gray', alpha=0.15, label=r'$\pm 2\sigma$ Poisson noise')
ax2.axhline(0, color='gray', linestyle=':')
ax2.set_xscale('log')
ax2.set_xlabel(r'Radius $r\; [\mathrm{kpc}]$')
ax2.set_ylabel(r'$(\rho - \bar{\rho}) / \bar{\rho}$')
ax2.set_ylim(-0.5, 0.5)
ax2.grid(True, which="both", ls=":", alpha=0.5)
ax2.legend(loc='lower left', fontsize=8)

plt.tight_layout()
plt.savefig('output/sanity_check_particles.png', dpi=300)
plt.show()