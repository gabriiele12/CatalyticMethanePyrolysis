"""Reaction-rate models for CH4 decomposition on Ni, Fe and activated carbon (AC).

Shared by reaction_rate.ipynb and Economics.ipynb.
Rates are in mol_CH4 / (kg_cat s); T is in K; pressures are given in bar and
converted to the unit each model uses.
"""
import numpy as np
from matplotlib.ticker import FuncFormatter, MultipleLocator

R = 8.314462618    # gas constant [J/(mol K)]
BAR_TO_PA = 1e5
ATM_TO_PA = 101325.0


# ---- Rate models ------------------------------------------------------------
def reaction_rate_nickel(T, P_CH4_bar, P_H2_bar):
    """Langmuir-Hinshelwood-type rate with CH4 and H2 adsorption (pressures in Pa)."""
    P_CH4 = P_CH4_bar * BAR_TO_PA
    P_H2 = P_H2_bar * BAR_TO_PA

    # Rate constant
    A, E = 9.77054, 88e3
    k = A * np.exp(-E / (R * T))

    # Adsorption constants
    A_CH4, deltaH_CH4 = 6.887737e-10, -56e3
    K_CH4 = A_CH4 * np.exp(-deltaH_CH4 / (R * T))
    A_H2, deltaH_H2 = 5.8009466e-16, -144e3
    K_H2 = A_H2 * np.exp(-deltaH_H2 / (R * T))

    # Equilibrium constant of CH4 <-> C + 2 H2
    Kp = 5.088e10 * np.exp(-91.2e3 / (R * T))

    numerator = k * (P_CH4 - P_H2**2 / Kp)
    denominator = (1 + K_CH4 * P_CH4 + K_H2 * P_H2**1.5) ** 2
    return numerator / denominator


def reaction_rate_iron(T, P_CH4_bar, P_H2_bar):
    """First order in CH4 with an approach-to-equilibrium term (pressures in atm)."""
    P_CH4 = P_CH4_bar * BAR_TO_PA / ATM_TO_PA
    P_H2 = P_H2_bar * BAR_TO_PA / ATM_TO_PA
    P0 = 1.0 * BAR_TO_PA / ATM_TO_PA  # standard pressure, 1 bar [atm]

    # Equilibrium constant of CH4 <-> C + 2 H2, and the approach to equilibrium
    deltaG = 89658.88 - 102.27 * T - 0.00428 * T**2 - 2499358.99 / T  # [J/mol]
    Keq = np.exp(-deltaG / (R * T))
    Qp = P_H2**2 / (P_CH4 * P0)
    eta = Qp / Keq

    # Rate constant (Arrhenius around a reference temperature)
    T_ref, kp_ref, Ea = 1073.0, 1.3, 160e3
    kp = kp_ref * np.exp(-(Ea / R) * (1.0 / T - 1.0 / T_ref))

    return kp * P_CH4 * (1.0 - eta)


def reaction_rate_carbon(T, P_CH4_bar):
    """Power law of order n in CH4 (pressure in atm); the original rate is
    divided by 60."""
    P_CH4 = P_CH4_bar * BAR_TO_PA / ATM_TO_PA

    n, kp0, Ea = 0.5, 1.4e10, 216e3
    kp = kp0 * np.exp(-Ea / (R * T))

    r_original = kp * P_CH4**n
    return r_original / 60.0


# ---- Rate vs. temperature ---------------------------------------------------
# Temperature range [°C] over which each model is used, and its plot color
temperature_ranges_C = {"Ni": (550, 700), "Fe": (700, 850), "AC": (800, 950)}
catalyst_colors = {"Ni": "#E69F00", "Fe": "#7B4AB5", "AC": "#D62728"}
model_names = {"Ni": "Ni model", "Fe": "Fe model", "AC": "AC model"}

temperature_label = r"Temperature / [$^\circ\mathrm{C}$]"
rate_label = (
    r"Reaction rate / "
    r"[$\mathrm{mol}_{\mathrm{CH_4}}\cdot"
    r"\mathrm{kg}_{\mathrm{cat}}^{\mbox{-}1}\cdot"
    r"\mathrm{s}^{\mbox{-}1}$]"
)


def rate_vs_temperature(catalyst, P_CH4_bar, P_H2_bar=0, n_points=300):
    """Temperatures [°C] over the catalyst's model range, and the rate there."""
    T_C = np.linspace(*temperature_ranges_C[catalyst], n_points)
    T_K = T_C + 273.15
    if catalyst == "Ni":
        rate = reaction_rate_nickel(T_K, P_CH4_bar, P_H2_bar)
    elif catalyst == "Fe":
        rate = reaction_rate_iron(T_K, P_CH4_bar, P_H2_bar)
    else:
        rate = reaction_rate_carbon(T_K, P_CH4_bar)
    return T_C, rate


def model_label(catalyst):
    """Legend label, e.g. 'Ni model'."""
    return model_names[catalyst]


def decimal_label(value, pos=None):
    """Label exact powers of ten as plain numbers (0.01, 0.1, 1, 10); others blank."""
    if value <= 0:
        return ""
    exponent = int(np.round(np.log10(value)))
    if not np.isclose(value, 10.0**exponent):
        return ""
    return f"{10.0**exponent:g}"


def draw_reaction_rates(ax, linewidth, P_CH4_bar=1.0, P_H2_bar=0):
    """Draw the rate of every catalyst vs. temperature on a log axis.

    Returns {catalyst: (T_C, rate)}. Axis labels, ticks and the legend are left
    to the caller, so each figure can use its own font sizes.
    """
    curves = {}
    for catalyst in temperature_ranges_C:
        T_C, rate = rate_vs_temperature(catalyst, P_CH4_bar, P_H2_bar)
        curves[catalyst] = (T_C, rate)
        # Only positive rates can be drawn on the log axis
        ax.plot(T_C, np.where(rate > 0, rate, np.nan), color=catalyst_colors[catalyst],
                linestyle="-", linewidth=linewidth, label=model_label(catalyst))

    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(decimal_label))

    # Temperature axis: 5 °C of padding on both sides so the curves don't touch the
    # frame; a labeled tick every 100 °C, minor ticks every 25 °C
    ax.set_xlim(545, 955)
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(25))
    return curves
