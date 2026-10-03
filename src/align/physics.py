"""Physical criteria shared by all phases (SI units)."""
import numpy as np

MU0 = 4e-7 * np.pi
KB = 1.380649e-23


def ms_ratio(T, Tc, s=0.5, p=2.5):
    """Kuz'min (PRL 94, 107204, 2005) shape of Ms(T)/Ms(0); s=0.5, p=5/2 used for every
    phase (the parameter only affects the strength estimate of a field bias, never a
    ceiling). Returns 0 above Tc."""
    t = np.clip(np.asarray(T, float) / Tc, 0, None)
    v = 1 - s * t ** 1.5 - (1 - s) * t ** p
    return np.where(t < 1, np.clip(v, 0, None) ** (1 / 3), 0.0)


def ms_at(T, Js300, Tc):
    """Ms (A/m) at T from the 300 K polarisation, via the same Ms(T) shape."""
    return Js300 / MU0 * ms_ratio(T, Tc) / ms_ratio(300.0, Tc)


def ferro_bias_density(K, Ms, H):
    """Energy density (J/m^3) favouring a variant with easy axis along H over one at 90 deg:
    K h (2 - h) for h = H/H_K < 1 and K for h >= 1 (Stoner-Wohlfarth minimum energies)."""
    HK = 2 * K / (MU0 * Ms) if Ms > 0 else np.inf
    h = H / HK
    return K * h * (2 - h) if h < 1 else K


def para_bias_density(chi, dchi_over_chi, H):
    """Paramagnetic anisotropy energy density mu0 * dchi * H^2 / 2."""
    return MU0 * chi * dchi_over_chi * H ** 2 / 2


def curie_weiss_chi(n_per_m3, mu_eff_muB, T, theta):
    """SI volume susceptibility C/(T - theta), C = mu0 n mu_eff^2 / (3 kB)."""
    muB = 9.2740100783e-24
    C = MU0 * n_per_m3 * (mu_eff_muB * muB) ** 2 / (3 * KB)
    return C / (T - theta)


def selection_fraction(x, n_variants=3):
    """Fraction of the favoured variant when the bias energy per deciding volume is x kT
    and the other n-1 variants are disfavoured by the full bias (Boltzmann)."""
    return np.exp(x) / (np.exp(x) + (n_variants - 1))
