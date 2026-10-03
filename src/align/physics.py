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


MU_B = 9.2740100783e-24


def mean_field_chi(Js0, n_mag_per_m3, Tc, T):
    """Mean-field estimate of the paramagnetic susceptibility above Tc from the
    saturation polarisation Js0 (T, ~0 K) and the density of magnetic atoms:
    moment m = Ms0/n (g=2, S=m/2), C = mu0 n g^2 S(S+1) muB^2 / (3 kB), chi = C/(T-Tc).
    An estimate, the same for every phase; used only to bound a paramagnetic bias."""
    m = Js0 / MU0 / n_mag_per_m3 / MU_B          # muB per magnetic atom
    S = m / 2
    C = MU0 * n_mag_per_m3 * 4 * S * (S + 1) * MU_B ** 2 / (3 * KB)
    return C / (T - Tc)


def para_bias_kT_bound(Js0, n_mag_per_m3, Tc, T, B_tesla, V_star):
    """Upper bound on the paramagnetic variant bias in units of kT:
    mu0 * chi * H^2 * V* / (2 kT) with dchi <= chi (uniaxial bound)."""
    chi = mean_field_chi(Js0, n_mag_per_m3, Tc, T)
    H = B_tesla / MU0
    return MU0 * chi * H ** 2 * V_star / (2 * KB * T)


def ferro_bias_kT(Js300, K300, Tc, T, B_tesla, V_star, n_callen=3):
    """Ferromagnetic variant bias (kT units) at T < Tc: Stoner-Wohlfarth energy
    difference between easy axis along and perpendicular to the field, with
    Ms(T) from the Kuz'min shape and K(T) = K300 (Ms(T)/Ms(300))^n (Callen-Callen n=3).
    Estimates; the ceiling only uses whether Tc > T."""
    if T >= Tc:
        return 0.0
    Ms = ms_at(T, Js300, Tc)
    K = K300 * (ms_ratio(T, Tc) / ms_ratio(300.0, Tc)) ** n_callen
    return ferro_bias_density(float(K), float(Ms), B_tesla / MU0) * V_star / (KB * T)
