import numpy as np
import matplotlib.pyplot as plt

def analytical_profile(r, rho_0=3.746e9, r_0=1.0):
    return rho_0 / (1.0 + (r / r_0)**2)

runs = {
    "NGP 100^3": ("output/profile_ngp_100.csv", "output/density_ngp_100.png"),
    "CIC 100^3": ("output/profile_cic_100.csv", "output/density_cic_100.png"),
    "NGP 200^3": ("output/profile_ngp_200.csv", "output/density_ngp_200.png"),
    "CIC 200^3": ("output/profile_cic_200.csv", "output/density_cic_200.png"),
}

r_fine = np.logspace(-1, np.log10(250), 300)
rho_exact = analytical_profile(r_fine)

for label, (csv_file, out_png) in runs.items():
    try:
        data = np.loadtxt(csv_file, delimiter=',', skiprows=1)
        r, rho = data[:, 0], data[:, 1]
    except FileNotFoundError:
        print(f"Skipping {csv_file} (not found)")
        continue

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True, gridspec_kw={'height_ratios': [2.5, 1]})

    # Top panel: densities
    ax1.plot(r_fine, rho_exact, 'k--', label='Analytical $\\rho(r)$', lw=1.5)
    ax1.plot(r, rho, label=label, color='crimson')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'$\rho(r)\; [M_\odot / \mathrm{kpc}^3]$')
    ax1.set_title(f'Density Profile: {label}')
    ax1.grid(True, which="both", ls=":", alpha=0.5)
    ax1.legend()

    # Bottom panel: fractional residuals
    rho_ref = analytical_profile(r)
    ax2.plot(r, (rho - rho_ref) / rho_ref, color='crimson')
    ax2.axhline(0, color='gray', linestyle=':')
    ax2.set_xscale('log')
    ax2.set_xlabel(r'Radius $r\; [\mathrm{kpc}]$')
    ax2.set_ylabel(r'$(\rho - \rho_{\mathrm{true}}) / \rho_{\mathrm{true}}$')
    ax2.set_ylim(-0.6, 0.6)
    ax2.grid(True, which="both", ls=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close(fig)  # Close to prevent memory accumulation in batch