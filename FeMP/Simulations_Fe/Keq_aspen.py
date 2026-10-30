import numpy as np

R_J = 8.3144598          # J/(mol K)
R_bar = 0.083144598      # bar m3/(kmol K)
P0 = 100000  #1 bar

def delta_g0(T):
    """
    Gibbs free energy variation [J/mol], T in K.
    Reaction: CH4 -> C(s) + 2H2
    """
    return 89658.88 - 102.27 * T - 0.00428 * (T ** 2) - 2499358.99 / T



def ln_inv_Kc_app(T):
    """
    Direct expression for ln(1/Kc_app).
    """
    n = 1

    return n* (delta_g0(T) / (R_J * T)) - n* (np.log(P0))

def ln_kp_aspen(T, A, B, C, D):
    return A + B / T + C * np.log(T) + D * T  #aspen equation for ln of Kp

def fit_kp_aspen():
    """
    Fit ln(Kp) to Aspen form:
        ln(Kp) = A + B/T + C*ln(T) + D*T

    Temperature range: 700–950 °C

    Returns
    -------
    A, B, C, D
    """
    # range temperature
    T_C = np.linspace(550, 1000, 451) #difference+1 to have 1C steps
    T_K = T_C + 273.15

    y = ln_inv_Kc_app(T_K)  #values of ln of Kp from paper definition

    X = np.column_stack([    #aspen coefficients according to its definition in description
        np.ones_like(T_K),  #coefficient A
        1.0 / T_K,  #coefficient B
        np.log(T_K),  #coefficietn C
        T_K  #coefficient D
    ])

    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)  #minimization of coefficients
    A, B, C, D = coeffs

    y_fit = ln_kp_aspen(T_K, A, B, C, D)   #values of ln of Kp from aspen definition

    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot

    return A, B, C, D