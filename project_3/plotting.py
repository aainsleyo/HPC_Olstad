# Packages imported
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
# Scipy with stat tests
from scipy.stats import kstest, shapiro, norm, gaussian_kde
import os
import warnings
# Tracy Widom distribution package
from TracyWidom import TracyWidom

warnings.filterwarnings('ignore')

# Tracy-Widom
tw      = TracyWidom(beta=1)
# set to one because that is the GOE setting
#now it will fit to curve
TW_X    = np.linspace(-6, 4, 1000)  
TW_PDF  = tw.pdf(TW_X)

#Paths
NARVAL_BASE = r"C:\Users\101040911\OneDrive - Ontario Tech University\Documents\HPC\project 3\results_narval"
WTABLE_PATH = (r"C:\Users\101040911\OneDrive - Ontario Tech University\Documents\HPC"
               r"\project 3\MCSC-6030G Project 3-20260414T201649Z-3-001"
               r"\MCSC-6030G Project 3\wtable")
EIGS_BASE   = (r"C:\Users\101040911\OneDrive - Ontario Tech University\Documents\HPC"
               r"\project 3\MCSC-6030G Project 3-20260414T201649Z-3-001"
               r"\MCSC-6030G Project 3")

# Loaders, Parse specific files to data files
def load_narval_eigs(n_sizes=(1000, 2000, 4000, 6000), ndat=1000):
    "large matrix sizes at 1000 iterations"
    data = {}
    for n in n_sizes:
        #in NARVAL_BASE path, find folder named "Eigs"
        f = os.path.join(NARVAL_BASE, f"n{n}_ndat{ndat}", "eigs")
        if not os.path.exists(f):
            print(f"  Missing: {f}"); continue
        raw = np.loadtxt(f)
        #only parsing the eigs, not the worker ID
        data[n] = raw[:, 1]
    return data

def load_connor_eigs(n_sizes=(1024, 2048, 4096, 8192), ndat=1024):
    "loading connor's data in to compare"
    data = {}
    for n in n_sizes:
        f = os.path.join(EIGS_BASE, f"eigs_{n}_{ndat}.out")
        if not os.path.exists(f):
            print(f"  Missing: {f}"); continue
        raw = np.loadtxt(f)
        data[n] = raw[:, 1]
    return data

# plotting the raw distributions
def plot_raw(data, fname, color='steelblue', title_suffix=''):
    n_sizes = sorted(data.keys())
    ncols = min(len(n_sizes), 2)
    nrows = (len(n_sizes) + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3.5*nrows))
    axes = np.array(axes).flatten()

    for ax, n in zip(axes, n_sizes):
        eigs = data[n]
        edge = 2 * np.sqrt(n)
        # the wigner semicircle law spectral edge for GOE matrices.
        # As n approaches infinity (large matrix size) the largest eigenvalue concentrates near 2√n.
        ax.hist(eigs, bins=40, density=True,
                color=color, edgecolor='white', alpha=0.85)
        ax.axvline(edge, color='crimson', linewidth=1.5,
                   label=r'$2\sqrt{n}$' + f' = {edge:.1f}')
        ax.set_xlim(np.percentile(eigs, 0.5), np.percentile(eigs, 99.9) * 1.005)
        #focus around the main histogram, this does not delete anything!
        ax.set_title(f'$n = {n}$', fontsize=11)
        ax.set_xlabel(r'$\lambda_{\max}$', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.legend(fontsize=9)

    for ax in axes[len(n_sizes):]:
        ax.set_visible(False)

    fig.suptitle(r'Empirical distribution of $\lambda_{\max}$' + title_suffix,
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(fname + '.pdf', bbox_inches='tight')
    plt.savefig(fname + '.png', dpi=150, bbox_inches='tight')
    plt.close()

# TW normalized
def plot_tw(data, fname, color='steelblue', title_suffix=''):
    n_sizes = sorted(data.keys())
    ncols = min(len(n_sizes), 2)
    nrows = (len(n_sizes) + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3.5*nrows))
    axes = np.array(axes).flatten()

    for ax, n in zip(axes, n_sizes):
        eigs = data[n]
        z    = (eigs - 2 * np.sqrt(n)) / (n ** (1/6))
        #TW normalization

        # display window only, no data removed
        z_display = z[(z > -0.6) & (z < 0.4)]
        ax.hist(z_display, bins=50, density=True,
                color=color, edgecolor='white', alpha=0.75, label='Empirical')

        kde = gaussian_kde(z_display)
        x   = np.linspace(-0.6, 0.4, 500)
        #window to see cluster clearer
        ax.plot(x, kde(x), color='royalblue' if color=='steelblue' else 'saddlebrown',
                linestyle=':', linewidth=1.5, label='KDE')

        # TW1 rescaled to empirical mean/std for shape comparison
        emp_mean = np.mean(z_display)
        emp_std  = np.std(z_display)
        TW1_MEAN, TW1_STD = -1.2065, np.sqrt(1.6078)
        x_tw_native = np.linspace(-6, 4, 500)
        x_tw_emp    = emp_mean + emp_std / TW1_STD * (x_tw_native - TW1_MEAN)
        tw_pdf_emp  = tw.pdf(x_tw_native) * TW1_STD / emp_std
        ax.plot(x_tw_emp, tw_pdf_emp, 'crimson', linewidth=2,
                label='TW$_1$ (shape)')

        # Gaussian fit
        mu, sigma = norm.fit(z_display)
        x_fit = np.linspace(-0.6, 0.4, 500)
        ax.plot(x_fit, norm.pdf(x_fit, mu, sigma),
                'forestgreen', linestyle='--', linewidth=1.5, label='Gaussian fit')

        ax.set_xlim(-0.6, 0.4)
        ax.set_title(f'$n = {n}$', fontsize=11)
        ax.set_xlabel(r'$z = (\lambda_{\max} - 2\sqrt{n})\,/\,n^{1/6}$', fontsize=9)
        ax.set_ylabel('Density', fontsize=10)
        ax.legend(fontsize=8)

    for ax in axes[len(n_sizes):]:
        ax.set_visible(False)

    fig.suptitle(r'Tracy–Widom normalized $\lambda_{\max}$' + title_suffix,
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(fname + '.pdf', bbox_inches='tight')
    plt.savefig(fname + '.png', dpi=150, bbox_inches='tight')
    plt.close()

# parallel scaling...
def plot_parallel_scaling(wtable_path):
    data   = np.loadtxt(wtable_path)
    n_vals = np.unique(data[:, 0]).astype(int)
    colors = plt.cm.tab10(np.linspace(0, 0.6, len(n_vals)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    for color, n in zip(colors, n_vals):
        mask  = data[:, 0] == n
        p_arr = data[mask, 1].astype(int)
        t_arr = data[mask, 2]
        idx   = np.argsort(p_arr)
        p_arr, t_arr = p_arr[idx], t_arr[idx]

        speedup    = t_arr[0] / t_arr
        efficiency = speedup / (p_arr / p_arr[0])

        ax1.plot(p_arr, speedup,          'o-', color=color, label=f'$n={n}$')
        ax2.plot(p_arr, efficiency * 100, 'o-', color=color, label=f'$n={n}$')

    p_range = np.array([data[:, 1].min(), data[:, 1].max()])
    ax1.plot(p_range, p_range / p_range[0], 'k--', linewidth=1.2, label='Ideal')
    ax2.axhline(100, color='k', linestyle='--', linewidth=1.2, label='Ideal (100%)')

    for ax, ylabel, title in zip(
        [ax1, ax2],
        ['Speedup $S(p)$', 'Efficiency $E(p)$ (%)'],
        ['Parallel Speedup', 'Parallel Efficiency']
    ):
        ax.set_xlabel('Number of workers $p$', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('plot_parallel_scaling.pdf', bbox_inches='tight')
    plt.savefig('plot_parallel_scaling.png', dpi=150, bbox_inches='tight')
    plt.close()

# Stats
def stats_table(data):
    hdr = f"  {'n':>6}  {'SW p':>9}  {'KS-N D':>8}  {'KS-N p':>9}  {'KS-TW D':>8}  {'Skew':>8}  {'Ex.Kurt':>8}"
    sep = "  " + "-" * (len(hdr) - 2)
    print(hdr)
    print(sep)

    for n in sorted(data.keys()):
        eigs = data[n]
        z    = (eigs - 2*np.sqrt(n)) / (n**(1/6))
        _, sw_p          = shapiro(z[:5000])
        ks_n, ks_n_p     = kstest((z - z.mean())/z.std(), 'norm')
        tw_cdf           = lambda x: tw.cdf(x)
        ks_tw, ks_tw_p   = kstest(z, tw_cdf)
        skew             = stats.skew(z)
        kurt             = stats.kurtosis(z)

        print(f"  {n:>6}  {sw_p:>9.2e}  {ks_n:>8.4f}  {ks_n_p:>9.2e}  {ks_tw:>8.4f}  {skew:>8.4f}  {kurt:>8.4f}")

    print(sep)

# Main to run all above
if __name__ == '__main__':
    narval = load_narval_eigs()
    if narval:
        plot_raw(narval, 'plot3_raw_narval')
        plot_tw(narval,  'plot4_tw_narval')
        stats_table(narval)

    narval_large = load_narval_eigs(n_sizes=(1000, 2000), ndat=10000)
    if narval_large:
        plot_raw(narval_large, 'plot3_raw_narval_large',
                 title_suffix=' ($M=10000$)')
        plot_tw(narval_large,  'plot4_tw_narval_large',
                title_suffix=' ($M=10000$)')
        stats_table(narval_large)

    connor = load_connor_eigs()
    if connor:
        plot_raw(connor, 'plot3_raw_connor', color='darkorange',
                 title_suffix=' (reference)')
        plot_tw(connor,  'plot4_tw_connor',  color='darkorange',
                title_suffix=' (reference)')
        stats_table(connor)

    if os.path.exists(WTABLE_PATH):
        plot_parallel_scaling(WTABLE_PATH)
    else:
        print(f"  wtable not found: {WTABLE_PATH}")
    print("Done.")
