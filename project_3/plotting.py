# Packages imported
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
# Scipy with stat tests
from scipy.stats import kstest, shapiro, norm, gaussian_kde
import os
import warnings

# Tracy Widom distribution package (important for GOE theory)
try:
    from TracyWidom import TracyWidom
except ImportError:
    print("ERROR: TracyWidom not installed. Run: pip install --user TracyWidom")
    exit(1)

warnings.filterwarnings('ignore')

# Tracy-Widom setup
tw = TracyWidom(beta=1)
# beta = 1 corresponds to GOE (real symmetric matrices)

# Paths SAFE for shell script
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# this is where all Narval output lives
NARVAL_BASE = os.path.join(BASE_DIR, "results_narval")

# parallel scaling table (if present)
WTABLE_PATH = os.path.join(BASE_DIR, "wtable")


# Loaders
def load_narval_eigs(n_sizes=(1000, 2000, 4000, 6000), ndat=1000):
    """
    Load eigenvalues from Narval runs.
    Each folder looks like: nXXXX_ndatYYYY/eigs
    """
    data = {}

    for n in n_sizes:
        # go into each folder and grab eigs file
        f = os.path.join(NARVAL_BASE, f"n{n}_ndat{ndat}", "eigs")

        if not os.path.exists(f):
            print(f"  Missing: {f}")
            continue

        raw = np.loadtxt(f)

        # only grabbing eigenvalues (column 1), not worker IDs
        data[n] = raw[:, 1]

    return data


# Raw distributions
def plot_raw(data, fname, color='steelblue', title_suffix=''):
    n_sizes = sorted(data.keys())

    ncols = min(len(n_sizes), 2)
    nrows = (len(n_sizes) + 1) // 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3.5*nrows))
    axes = np.array(axes).flatten()

    for ax, n in zip(axes, n_sizes):
        eigs = data[n]

        # Wigner semicircle edge for GOE
        # largest eigenvalue should sit near this
        edge = 2 * np.sqrt(n)

        ax.hist(eigs, bins=40, density=True,
                color=color, edgecolor='white', alpha=0.85)

        # vertical line showing theoretical edge
        ax.axvline(edge, color='crimson', linewidth=1.5,
                   label=r'$2\sqrt{n}$')

        # zoom in around where data actually is
        ax.set_xlim(np.percentile(eigs, 0.5),
                    np.percentile(eigs, 99.9) * 1.005)

        ax.set_title(f'$n = {n}$')
        ax.set_xlabel(r'$\lambda_{\max}$')
        ax.set_ylabel('Density')
        ax.legend()

    # hide unused axes
    for ax in axes[len(n_sizes):]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig(fname + '.pdf')
    plt.savefig(fname + '.png', dpi=150)
    plt.close()


# TW normalized plots
def plot_tw(data, fname, color='steelblue', title_suffix=''):
    n_sizes = sorted(data.keys())

    ncols = min(len(n_sizes), 2)
    nrows = (len(n_sizes) + 1) // 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3.5*nrows))
    axes = np.array(axes).flatten()

    for ax, n in zip(axes, n_sizes):
        eigs = data[n]

        # THIS is the Tracy–Widom scaling
        # subtract the edge and rescale fluctuations
        z = (eigs - 2 * np.sqrt(n)) / (n ** (1/6))

        # just zooming in visually, not removing data from stats
        z_display = z[(z > -0.6) & (z < 0.4)]

        ax.hist(z_display, bins=50, density=True,
                color=color, edgecolor='white', alpha=0.75,
                label='Empirical')

        # KDE just helps visualize shape
        kde = gaussian_kde(z_display)
        x = np.linspace(-0.6, 0.4, 500)
        ax.plot(x, kde(x), linestyle=':', linewidth=1.5,
                label='KDE')

        # Tracy–Widom theoretical curve
        # rescaled to match empirical mean/std to compare SHAPE only
        emp_mean = np.mean(z_display)
        emp_std  = np.std(z_display)

        TW1_MEAN, TW1_STD = -1.2065, np.sqrt(1.6078)

        x_tw_native = np.linspace(-6, 4, 500)

        x_tw_emp = emp_mean + emp_std / TW1_STD * (x_tw_native - TW1_MEAN)
        tw_pdf_emp = tw.pdf(x_tw_native) * TW1_STD / emp_std

        ax.plot(x_tw_emp, tw_pdf_emp, 'crimson', linewidth=2,
                label='TW$_1$ (shape)')

        # Gaussian fit
        mu, sigma = norm.fit(z_display)
        ax.plot(x, norm.pdf(x, mu, sigma),
                'forestgreen', linestyle='--', linewidth=1.5,
                label='Gaussian')

        ax.set_xlim(-0.6, 0.4)
        ax.set_title(f'$n = {n}$')
        ax.set_xlabel(r'$z = (\lambda_{\max} - 2\sqrt{n})/n^{1/6}$')
        ax.set_ylabel('Density')
        ax.legend()

    for ax in axes[len(n_sizes):]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig(fname + '.pdf')
    plt.savefig(fname + '.png', dpi=150)
    plt.close()


# Stats Table
def stats_table(data):
    print(" n     SW p      KS-N D   KS-N p    KS-TW D   Skew     Kurt")

    for n in sorted(data.keys()):
        eigs = data[n]

        # same TW normalization as above
        z = (eigs - 2*np.sqrt(n)) / (n**(1/6))

        # normality test
        _, sw_p = shapiro(z[:5000])

        # KS vs normal
        ks_n, ks_n_p = kstest((z - z.mean())/z.std(), 'norm')

        # KS vs Tracy–Widom
        ks_tw, _ = kstest(z, lambda x: tw.cdf(x))

        skew = stats.skew(z)
        kurt = stats.kurtosis(z)

        print(f"{n:4d}  {sw_p:9.2e}  {ks_n:8.4f}  {ks_n_p:9.2e}  "
              f"{ks_tw:8.4f}  {skew:8.4f}  {kurt:8.4f}")


# Main
if __name__ == '__main__':

    # make sure data exists
    if not os.path.exists(NARVAL_BASE):
        print("ERROR: results_narval folder not found.")
        exit(1)

    # standard runs (M = 1000)
    narval = load_narval_eigs()

    if narval:
        plot_raw(narval, 'plot_raw')
        plot_tw(narval,  'plot_tw')
        stats_table(narval)

    # larger ensemble (better tails)
    narval_large = load_narval_eigs(n_sizes=(1000, 2000), ndat=10000)

    if narval_large:
        plot_raw(narval_large, 'plot_raw_large',
                 title_suffix=' ($M=10000$)')
        plot_tw(narval_large,  'plot_tw_large',
                title_suffix=' ($M=10000$)')
        stats_table(narval_large)

    # parallel scaling (if file exists)
    if os.path.exists(WTABLE_PATH):
        data = np.loadtxt(WTABLE_PATH)
        print("Loaded parallel scaling data.")
    else:
        print(f"wtable not found: {WTABLE_PATH}")

    print("Done.")
