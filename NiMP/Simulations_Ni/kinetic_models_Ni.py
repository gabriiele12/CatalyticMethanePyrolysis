import time
import warnings
import math
import numpy as np
from scipy.optimize import minimize, root_scalar


def Arrhenius(A, E, T):
    """
    Arrhenius equation to compute temperature-dependent rate constants.

    Parameters:
        A (float): Pre-exponential factor.
        E (float): Activation energy (kJ/mol).
        T (float or ndarray): Temperature in Kelvin.

    Returns:
        float or ndarray: Rate constant(s).
    """
    R = 8.3144598e-3  # Universal gas constant in kJ / (mol·K)
    return A * np.exp(-E / (R * T))


def reaction_adsorption_constant(T, r_constants):
    """
    Calculate reaction and adsorption constants from temperature.

    Parameters:
        T (float): Temperature in Celsius.
        r_constants (dict): Contains pre-exponential factors and activation energies for:
            - 'k': reaction rate constant
            - 'K_p': equilibrium constant
            - 'K_ch4': CH₄ adsorption constant
            - 'K_h2': H₂ adsorption constant

    Returns:
        tuple: (k, Kp, k_ch4, k_h2)
    """
    T_k = T + 273.15  # Convert to Kelvin
    arrhenius = lambda A, E: Arrhenius(A, E, T_k)

    k = arrhenius(*r_constants['k'])
    Kp = arrhenius(*r_constants['K_p'])
    k_ch4 = arrhenius(*r_constants['K_ch4'])
    k_h2 = arrhenius(*r_constants['K_h2'])

    return k, Kp, k_ch4, k_h2


def deactivation_constant(T, da_constants):
    """
    Calculate deactivation-related constants based on temperature.

    Parameters:
        T (float): Temperature in Celsius.
        da_constants (dict): Contains Arrhenius parameters for deactivation mechanisms:
            - 'kd', 'kd_c', 'kd_ch4', 'kd_h2'

    Returns:
        tuple: (kd, kd_c, kd_ch4, kd_h2)
    """
    T += 273.15  # Convert to Kelvin
    kd = Arrhenius(*da_constants['kd'], T)
    kd_c = Arrhenius(*da_constants['kd_c'], T)
    kd_ch4 = Arrhenius(*da_constants['kd_ch4'], T)
    kd_h2 = Arrhenius(*da_constants['kd_h2'], T)

    return kd, kd_c, kd_ch4, kd_h2


def residenceTime(cat_holdup, cat_rate):
    """
    Compute residence time from catalyst holdup and catalyst rate.

    Parameters:
        cat_holdup (float): Catalyst weight (kg).
        cat_rate (float): Catalyst rate (kg/s).

    Returns:
        float: Residence time (s).
    """
    return cat_holdup / cat_rate


def deactivation(T, Pch4, Ph2, t, da_constants):
    """
    Evaluate catalyst deactivation factor over time.

    Parameters:
        T (float): Temperature in Celsius.
        Pch4 (float): Partial pressure of CH₄ (Pa).
        Ph2 (float): Partial pressure of H₂ (Pa).
        t (float or ndarray): Time(s) in seconds.
        da_constants (dict): Deactivation constant parameters.

    Returns:
        float or ndarray: Deactivation factor(s) (0 to 1).
    """
    kd, kd_c, kd_ch4, kd_h2 = deactivation_constant(T, da_constants)
    t = np.array(t)

    # Compute deactivation denominator safely
    denominator = 1 - (0.5 * kd * (kd_c + kd_ch4 * Pch4 + kd_h2 * (Ph2 ** 0.83)) * t)
    a = np.zeros_like(denominator)

    valid = denominator > 0
    a[valid] = (1 / denominator[valid]) ** (-0.8)

    a = np.clip(a, 0, 1)
    return a if t.size > 1 else a.item()


def RK4constants(T, Pch4, Ph2, prev_a, da_constants, dt):
    """
    Compute Runge-Kutta (RK4) constants for deactivation factor.

    Parameters:
        T (float): Temperature in Celsius.
        Pch4 (float): Partial pressure of CH₄ (Pa).
        Ph2 (float): Partial pressure of H₂ (Pa).
        prev_a (float): Previous deactivation factor.
        da_constants (dict): Deactivation constant parameters.
        dt (float): Time step (s).

    Returns:
        tuple: RK4 constants (k1, k2, k3, k4) for use in integration.
    """
    kd, kd_c, kd_ch4, kd_h2 = deactivation_constant(T, da_constants)

    # Estimate time where deactivation matches previous a
    def equation_to_solve(t, Pch4, Ph2, kd, kd_c, kd_ch4, kd_h2, prev_a):
        term = min(0.5 * kd * t * (kd_c + kd_ch4 * Pch4 + kd_h2 * (Ph2 ** 0.83)), 1)
        return min((1 - term) ** 0.8, 1) - prev_a

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Derivative was zero")
        sol = root_scalar(equation_to_solve, args=(Pch4, Ph2, kd, kd_c, kd_ch4, kd_h2, prev_a), x0=0, method='newton')

    t0 = sol.root if sol.converged else 0

    # Derivative of deactivation with respect to time
    def dadt(t, Pch4, Ph2, kd, kd_c, kd_ch4, kd_h2):
        factor = Pch4 * kd_ch4 + (Ph2 ** 0.83) * kd_h2 + kd_c
        numerator = -0.4 * kd * factor
        term = min(0.5 * kd * t * factor, 1)
        denominator = ((1 - term) ** 0.2) + 1e-12
        return min(numerator / denominator, 0)

    # Compute RK4 terms
    k1 = dt * dadt(t0, Pch4, Ph2, kd, kd_c, kd_ch4, kd_h2)
    k2 = dt * dadt(t0 + dt / 2, Pch4, Ph2, kd, kd_c, kd_ch4, kd_h2)
    k3 = dt * dadt(t0 + dt / 2, Pch4, Ph2, kd, kd_c, kd_ch4, kd_h2)
    k4 = dt * dadt(t0 + dt, Pch4, Ph2, kd, kd_c, kd_ch4, kd_h2)

    return k1, k2, k3, k4


def initial_reaction_rate(T, Pch4, Ph2, r_constants):
    """
    Compute initial reaction rate based on Langmuir-Hinshelwood kinetics.

    Parameters:
        T (float): Temperature in Celsius.
        Pch4 (float): Partial pressure of CH₄ (Pa).
        Ph2 (float): Partial pressure of H₂ (Pa).
        r_constants (dict): Reaction and adsorption constants.

    Returns:
        float: Initial reaction rate (mol/kg·s).
    """
    k, Kp, k_ch4, k_h2 = reaction_adsorption_constant(T, r_constants)
    adsorption = (1 + k_h2 * (Ph2 ** 1.5) + k_ch4 * Pch4) ** 2
    driving_force = Pch4 - (Ph2 ** 2) / Kp
    r0 = k * driving_force / adsorption

    return r0



def isothermal_cstr(sim_config, Fch4_in, Fh2_in, p_drop, cat_wt, init_guess, verbose=False):
    """
    Simulate the isothermal CMP-FBR Continuous Stirred Tank Reactor (CSTR) model.

    The reactor is evaluated under two operation modes:
    - 'constant_rate': maintains a fixed catalyst rate.
    - 'constant_da': adjusts catalyst rate to achieve a target average deactivation factor.

    Parameters:
        sim_config (dict): Reactor configuration, including temperature, pressure, kinetics, deactivation, and mode
        settings.
        Fch4_in (float): Inlet molar flow rate of CH₄ (mol/s).
        Fh2_in (float): Inlet molar flow rate of H₂ (mol/s).
        p_drop (float): Pressure drop across the reactor (bar).
        cat_wt (float): Catalyst weight (kg).
        init_guess (float): Initial guess for catalyst rate (used in 'constant_da' mode).
        verbose (bool): If True, displays simulation details.

    Returns:
        dict: Simulation results, including outlet flow rates, pressures, average deactivation, catalyst rate, and
        residence time.
    """
    Pch4, Ph2, Fch4_out, Fh2_out, Fc_out, average_a, cat_rate, tres = [None] * 8
    # Extract key configuration parameters
    temp = sim_config['reactor']['temperature']             # Reactor temperature (°C)
    press = (sim_config['comp1_pressure']-p_drop) * 1e5     # Convert pressure (bar to Pa)
    mode_type = sim_config['mode']['type']                  # Operating mode type
    mode_value = sim_config['mode']['value']                # Target value for mode
    r_constants = sim_config['reaction_constants']          # Reaction rate constants
    da_constants = sim_config['deactivation_constants']     # Catalyst deactivation constants
    step_size = sim_config['mode']['step_size']             # Adjustment step size (for constant_da mode)

    def log(message):
        """Print progress message if verbose mode is enabled."""
        if verbose:
            print(message)

    def F(t, tau):
        """
        Cumulative distribution function for RTD of ideal CSTR.

        Parameters:
            t (float): Time (s).
            tau (float): Mean residence time (s).

        Returns:
            float: Cumulative fraction exited.
        """
        return 1 - math.exp(-t / tau)

    def RTD(tau):
        """
        Cumulative distribution function for RTD of ideal CSTR.

        Parameters:
            tau (float): Mean residence time (s).

        Returns:
            float: Cumulative fraction exited.
        """
        fractions = []
        times = []
        dt = 0

        while dt < (7 * tau):
            result = F(dt + (tau / 10), tau) - F(dt, tau)
            fractions.append(result)
            times.append(dt + (tau / 10) - (tau / 20))
            dt += tau / 10

        fractions = np.array(fractions)
        fractions /= np.sum(fractions)
        return fractions, np.array(times)

    def reactor_equations(flows, Fch4_in, Fh2_in, temp, press, r_constants, da_constants, cat_wt, tres):
        """
        Mass balance equations for CH₄ and H₂ in a CSTR.

        Returns:
            float: Sum of squared residuals for CH₄ and H₂ balances.
            float: Average deactivation factor.
        """
        Fch4_out, Fh2_out = flows
        Ftot = Fch4_out + Fh2_out
        Pch4 = Fch4_out * press / Ftot
        Ph2 = press - Pch4

        r0 = initial_reaction_rate(temp, Pch4, Ph2, r_constants)
        particle_fractions, time_in_reactor = RTD(tres)
        average_a = np.sum(deactivation(temp, Pch4, Ph2, time_in_reactor, da_constants)*particle_fractions)
        r = r0 * average_a

        R_ch4 = Fch4_in - Fch4_out - r * cat_wt
        R_h2 = Fh2_in - Fh2_out + 2 * r * cat_wt

        return R_ch4 ** 2 + R_h2 ** 2, average_a

    def simulate(cat_rate, Fch4_in, Fh2_in, temp, press, r_constants, da_constants):
        """
        Solve the CSTR model for a given catalyst rate.

        Returns:
            tuple: Outlet flow rates, partial pressures, deactivation, catalyst rate, and residence time.
        """
        tres = residenceTime(cat_wt, cat_rate)
        init_value = np.array([50, 10])  # Initial guess for outlet flow rates

        # Minimize mass balance residuals to find outlet flows
        result = minimize(lambda flows:
                          reactor_equations(flows, Fch4_in, Fh2_in, temp, press,
                                            r_constants, da_constants, cat_wt, tres)[0],
                          x0=init_value,
                          bounds=[(1e-11, 1e5), (1e-11, 1e5)])

        Fch4_out, Fh2_out = result.x
        Ftot = Fch4_out + Fh2_out
        Pch4 = Fch4_out * press / Ftot
        Ph2 = press - Pch4

        particle_fractions, time_in_reactor = RTD(tres)
        average_a = float(np.sum(deactivation(temp, Pch4, Ph2, time_in_reactor, da_constants)*particle_fractions))

        if verbose:
            log(f"Simulation complete: a = {average_a:.4f}, Pch4 = {Pch4 / 1e5:.3f} bar, Ph2 = {Ph2 / 1e5:.3f} bar")

        Fc_out = initial_reaction_rate(temp, Pch4, Ph2, r_constants) * average_a * cat_wt

        return Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres

    # Mode: constant catalyst rate
    if mode_type == 'constant_rate':
        cat_rate = mode_value
        Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres = simulate(cat_rate, Fch4_in, Fh2_in, temp,
                                                                                   press, r_constants, da_constants)

    # Mode: iterate catalyst rate to achieve target average deactivation
    elif mode_type == 'constant_da':
        error = 1
        cat_rate = init_guess

        while error > 0.001:
            Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres = simulate(cat_rate, Fch4_in, Fh2_in, temp,
                                                                                       press, r_constants, da_constants)
            print(f"Running isothermal CSTR model. DA:{average_a:.4g}; CatRate:{cat_rate:.4g}", end='\r', flush=True)
            time.sleep(0.1)

            # Adjust catalyst rate based on error direction
            error = abs(average_a - mode_value)
            if average_a > mode_value:
                cat_rate -= error * init_guess * step_size
            elif average_a < mode_value:
                cat_rate += error * init_guess * step_size

            if cat_rate < 0:
                raise Exception("Catalyst rate has reached 0")

    # Return simulation results
    return {
        "Pch4": Pch4 / 1e5,  # Convert to bar
        "Ph2": Ph2 / 1e5,  # Convert to bar
        "Fch4_out": Fch4_out,
        "Fh2_out": Fh2_out,
        "Fc_out": Fc_out,
        "average_a": average_a,
        "cat_rate": cat_rate,
        "tres": tres
    }


def isothermal_pfr(sim_config, Fch4_in, Fh2_in, p_drop, cat_wt, init_guess, n=501, verbose=False):
    """
    Simulate the isothermal CMP-FBR Plug Flow Reactor (PFR) model.

    The reactor is evaluated under two operation modes:
    - 'constant_rate': maintains a fixed catalyst rate.
    - 'constant_da': adjusts catalyst rate to achieve a target average deactivation factor.

    Parameters:
        sim_config (dict): Configuration dictionary with reactor conditions, kinetic and deactivation parameters.
        Fch4_in (float): Inlet molar flow rate of CH₄ (mol/s).
        Fh2_in (float): Inlet molar flow rate of H₂ (mol/s).
        p_drop (float): Total pressure drop across reactor (bar).
        cat_wt (float): Catalyst weight (kg).
        init_guess (float): Initial guess for catalyst rate (used in constant_da mode).
        n (int): Number of discretization segments along reactor length.
        verbose (bool): If True, print simulation progress during iteration.

    Returns:
        dict: Simulation results, including outlet flow rates, pressures, average deactivation, catalyst rate, and
        residence time.
    """

    # Initialize result placeholders
    Pch4, Ph2, Fch4_out, Fh2_out, Fc_out, average_a, cat_rate, tres = [None] * 8

    # Extract parameters from configuration
    temp = sim_config['reactor']['temperature']            # Reactor temperature (°C)
    press = sim_config['comp1_pressure'] * 1e5             # Inlet pressure (bar to Pa)
    mode_type = sim_config['mode']['type']                 # Operating mode type
    mode_value = sim_config['mode']['value']               # Target value for mode
    r_constants = sim_config['reaction_constants']         # Reaction kinetic constants
    da_constants = sim_config['deactivation_constants']    # Catalyst deactivation parameters
    step_size = sim_config['mode']['step_size']            # Tuning step size (constant_da mode)

    def log(message):
        """Print progress message if verbose is enabled."""
        if verbose:
            print(message)

    def update_deactivation(temp, Pch4, Ph2, a_prev, da_constants, dt, initial=False):
        """
        Compute deactivation factor at a segment using either initial deactivation expression or RK4 integration.

        Parameters:
            temp (float): Reactor temperature (°C).
            Pch4 (float): CH₄ partial pressure (Pa).
            Ph2 (float): H₂ partial pressure (Pa).
            a_prev (float): Previous segment's deactivation factor.
            da_constants (dict): Deactivation parameters.
            dt (float): Time increment per reactor segment (s).
            initial (bool): Whether to use the initial deactivation formula.

        Returns:
            float: Updated deactivation factor (bounded between 0 and 1).
        """
        if initial:
            a = deactivation(temp, Pch4, Ph2, dt, da_constants)
        else:
            k1, k2, k3, k4 = RK4constants(temp, Pch4, Ph2, a_prev, da_constants, dt)
            a = a_prev + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        return max(min(a, 1), 0)

    def reactor_equations(flows, Fch4_in, Fh2_in, temp, press, r_constants, da_constants, cat_wt, n, dt, a_prev, i):
        """
        Define the objective function for optimization at each reactor segment.

        Parameters:
            flows (array): Current CH₄ and H₂ flow rates to be optimized.
            (Other parameters are used to calculate the reaction rate and update deactivation.)

        Returns:
            float: Squared residuals of CH₄ and H₂ balances.
        """
        Fch4_out, Fh2_out = flows
        Ftot = Fch4_out + Fh2_out
        Pch4 = Fch4_out * press / Ftot
        Ph2 = press - Pch4

        a = update_deactivation(temp, Pch4, Ph2, a_prev, da_constants, dt, initial=(i == 0))

        r = initial_reaction_rate(temp, Pch4, Ph2, r_constants) * a

        # Residuals from molar balances
        R_ch4 = Fch4_in - Fch4_out - r * cat_wt / n
        R_h2 = Fh2_in - Fh2_out + 2 * r * cat_wt / n

        return R_ch4 ** 2 + R_h2 ** 2

    def simulate(cat_rate, Fch4_in, Fh2_in, temp, press, p_drop, r_constants, da_constants):
        """
        Perform simulation over the full reactor length with given catalyst rate.

        Returns:
            tuple: Outlet flows, pressures, average deactivation, and other results.
        """
        Fch4_out, Fh2_out, Fc_out, Pch4, Ph2 = [None] * 5
        Fc_in = 0
        tres = residenceTime(cat_wt, cat_rate)
        dt = tres / n
        a = np.zeros(n)
        a_prev = 1
        init_value = np.array([50, 10])  # Initial guess for flow optimization

        for i in range(n):
            result = minimize(fun=reactor_equations,
                              x0=init_value,
                              args=(Fch4_in, Fh2_in, temp, press, r_constants, da_constants, cat_wt, n, dt, a_prev, i),
                              bounds=[(1e-11, 100000), (1e-11, 100000)])
            Fch4_out, Fh2_out = result.x
            Ftot = Fch4_out + Fh2_out
            Pch4 = Fch4_out * press / Ftot
            Ph2 = press - Pch4

            a[i] = update_deactivation(temp, Pch4, Ph2, a_prev, da_constants, dt, initial=(i == 0))
            a_prev = a[i]
            r = initial_reaction_rate(temp, Pch4, Ph2, r_constants) * a[i]
            Fc_out = Fc_in + r * cat_wt / n

            # Update inlet flows and pressure for the next segment
            Fch4_in, Fh2_in, Fc_in = Fch4_out, Fh2_out, Fc_out
            press -= p_drop * 1e5 / n

            if verbose and i % 100 == 0:
                log(f"[Step {i}/{n}] a = {a[i]:.4f}, Pch4 = {Pch4 / 1e5:.3f} bar, Ph2 = {Ph2 / 1e5:.3f} bar")

        return Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, float(np.average(a)), tres

    # Mode: constant catalyst rate
    if mode_type == 'constant_rate':
        cat_rate = mode_value
        Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres = simulate(cat_rate, Fch4_in, Fh2_in, temp,
                                                                                   press, p_drop, r_constants,
                                                                                   da_constants)

    # Mode: adjust catalyst rate to meet a target deactivation factor
    elif mode_type == 'constant_da':
        error = 1
        cat_rate = init_guess
        max_iter = 700
        iteration = 0

        while error > 0.001 and iteration < max_iter:

            Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres = simulate(cat_rate, Fch4_in, Fh2_in, temp,
                                                                                       press, p_drop, r_constants,
                                                                                       da_constants)
            print(f"Running isothermal PFR model. DA:{average_a:.4g}; CatRate:{cat_rate:.4g}", end='\r', flush=True)

            # Adjust catalyst rate based on error direction
            error = abs(average_a - mode_value)
            if average_a > mode_value:
                cat_rate -= error * init_guess * step_size
            elif average_a < mode_value:
                cat_rate += error * init_guess * step_size

            if cat_rate < 0:
                raise Exception("Catalyst rate has reached 0")

            iteration += 1

        if iteration == max_iter:
            raise RuntimeError("PFR constant_da did not converge")

    # Return simulation results
    return {
        "Pch4": Pch4 / 1e5,  # Convert to bar
        "Ph2": Ph2 / 1e5,  # Convert to bar
        "Fch4_out": Fch4_out,
        "Fh2_out": Fh2_out,
        "Fc_out": Fc_out,
        "average_a": average_a,
        "cat_rate": cat_rate,
        "tres": tres
    }