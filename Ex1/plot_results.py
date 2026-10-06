import numpy as np
import matplotlib.pyplot as plt

def analytical_profile(r, rho_0=3.746e9, r_0=1.0):
    return rho_0 / (1.0 + (r / r_0)**2)

# Load CSV files exported by C++
runs = {
    "NGP 100^3": "output/profile_ngp_100.csv",
    "CIC 100^3": "output/profile_cic_100.csv",
    "NGP 200^3": "output/profile_ngp_200.csv",
    "CIC 200^3": "output/profile_cic_200.csv",
}

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9), sharex=True, gridspec_kw={'height_ratios': [2.5, 1]})

r_fine = np.logspace(-1, np.log10(250), 300)
rho_exact = analytical_profile(r_fine)
ax1.plot(r_fine, rho_exact, 'k--', label='Analytical $\\rho(r)$', lw=2)

for label, filename in runs.items():
    try:
        data = np.loadtxt(filename, delimiter=',', skiprows=1)
        r, rho = data[:, 0], data[:, 1]
        
        # Density profile
        line, = ax1.plot(r, rho, label=label, alpha=0.85)
        
        # Fractional error relative to exact profile
        rho_ref = analytical_profile(r)
        ax2.plot(r, (rho - rho_ref) / rho_ref, label=label, color=line.get_color(), alpha=0.85)
    except FileNotFoundError:
        print(f"Skipping {filename} (not yet generated)")

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$\rho(r)\; [M_\odot / \mathrm{kpc}^3]$')
ax1.set_title('Galaxy Density Profile: NGP vs CIC')
ax1.grid(True, which="both", ls=":", alpha=0.5)
ax1.legend()

ax2.set_xscale('log')
ax2.set_xlabel(r'Radius $r\; [\mathrm{kpc}]$')
ax2.set_ylabel(r'$(\rho_{\mathrm{grid}} - \rho_{\mathrm{true}}) / \rho_{\mathrm{true}}$')
ax2.axhline(0, color='gray', linestyle=':')
ax2.set_ylim(-0.6, 0.6)
ax2.grid(True, which="both", ls=":", alpha=0.5)

plt.tight_layout()
plt.savefig('output/density_comparison.png', dpi=300)
plt.show()