import warnings
import math
import numpy as np
from scipy.optimize import minimize, root_scalar


def Arrhenius(A, E, T, T_reference):
    """
    Arrhenius equation to compute temperature-dependent rate constants.

    Parameters:
        A (float): Pre-exponential factor.
        E (float): Activation energy (kJ/mol).
        T (float or ndarray): Temperature in Kelvin.

    Returns:
        float or ndarray: Rate constant(s).
    """
    R = 8.3144598e-3 # Universal gas constant in kJ / (mol·K)
    return A * np.exp((-E / R) * (1/T - 1/T_reference))


def reaction_adsorption_constant(T, r_constants):
    """
    Calculate the temperature-dependent reaction rate constant and equilibrium constant.

    Parameters:
        T (float): Temperature in Celsius.
        r_constants (dict): Contains the pre-exponential factor, activation energy, and reference temperature for 'k'.

    Returns:
        tuple: (k, Keq), where k is the reaction rate constant and Keq is the dimensionless equilibrium constant.
    """
    T_k = T + 273.15  # Convert to Kelvin
    A, E, T_ref = r_constants["k"]
    k = Arrhenius(A, E, T_k, T_ref)

    delta_G0 = 89658.88 - 102.27 * T_k - 0.00428 * (T_k ** 2) - 2499358.99 / T_k  #in J/mol
    R = 8.3144598  # Universal gas constant in J / (mol·K)
    Keq = np.exp(-delta_G0 / (R * T_k))
    return k, Keq


def deactivation_constant(T, da_constants):
    T_k = T + 273.15

    A_d, E_d, T_refd = da_constants["kd"]

    kd = Arrhenius(A_d, E_d, T_k, T_refd)
    return kd


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



def deactivation(T, Pch4, Ph2, t, r_constants, da_constants):
    kd = deactivation_constant(T, da_constants)
    r0 = initial_reaction_rate(T, Pch4, Ph2, r_constants)

    a_inf = 0.001
    M_C = 12.011  # kg/kmol

    t = np.asarray(t, dtype=float)

    coefficient = kd * M_C * r0

    denominator = (1 - (1 - a_inf) * np.exp(-a_inf * coefficient * t))

    a = a_inf / denominator
    a = np.clip(a, a_inf, 1.0)

    return a.item() if a.ndim == 0 else a


def RK4constants(T, Pch4, Ph2, a_prev, r_constants, da_constants, dt):
    kd = deactivation_constant(T, da_constants)

    a_inf = 0.001
    M_C = 12.011  # kgc/kmolc = g/mol

    r0 = initial_reaction_rate(T, Pch4, Ph2, r_constants)
    coefficient = kd * M_C * r0

    # Exact solution by inverting the equation
    t0 = np.log(a_prev * (1 - a_inf) / (a_prev - a_inf)) / (a_inf * coefficient)

    # Derivative with respect to time
    def dadt(t, coefficient):
        denominator = 1 - (1 - a_inf) * np.exp(-a_inf * coefficient * t)
        a_t = a_inf / denominator
        derivative = -coefficient * a_t * (a_t - a_inf)
        return derivative

    #RK4
    k1 = dt * dadt(t0, coefficient)
    k2 = dt * dadt(t0 + dt / 2, coefficient)
    k3 = dt * dadt(t0 + dt / 2, coefficient)
    k4 = dt * dadt(t0 + dt, coefficient)

    return k1, k2, k3, k4


def initial_reaction_rate(T, Pch4, Ph2, r_constants):
    """
    Compute the initial methane-pyrolysis rate using the reversible first-order kinetic expression.

    Parameters:
        T (float): Temperature in Celsius.
        Pch4 (float): Partial pressure of CH₄ (Pa).
        Ph2 (float): Partial pressure of H₂ (Pa).
        r_constants (dict): Reaction-rate parameters used to calculate k and Keq.

    Returns:
        float: Initial reaction rate (kmol/(kg_cat·s)).
    """
    k, Keq = reaction_adsorption_constant(T, r_constants)
    n = 2
    P0 = 100000 # 1 bar
    r0 = k * (Pch4 - (Ph2 ** n)/(Keq * P0))

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
    press = (sim_config['comp1_pressure']-p_drop) * 1e5     # Convert pressure from bar to Pa (simulation and Aspen work in bar, here in Pacal)
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
        average_a = np.sum(deactivation(temp, Pch4, Ph2, time_in_reactor, r_constants, da_constants)*particle_fractions)
        r = r0 * average_a

        R_ch4 = Fch4_in - Fch4_out - r * cat_wt * 1000
        R_h2 = Fh2_in - Fh2_out + 2 * r * cat_wt * 1000


        return R_ch4 ** 2 + R_h2 ** 2, average_a

    def simulate(cat_rate, Fch4_in, Fh2_in, temp, press, r_constants, da_constants):
        """
        Solve the CSTR model for a given catalyst rate.

        Returns:
            tuple: Outlet flow rates, partial pressures, deactivation, catalyst rate, and residence time.
        """
        tres = residenceTime(cat_wt, cat_rate)
        init_value = np.array([500, 100])  # Initial guess for outlet flow rates

        # Minimize mass balance residuals to find outlet flows
        result = minimize(lambda flows:
                          reactor_equations(flows, Fch4_in, Fh2_in, temp, press,
                                            r_constants, da_constants, cat_wt, tres)[0],
                          x0=init_value,
                          bounds=[(1e-11, 1e6), (1e-11, 1e6)])

        if not result.success or result.fun > 1e-8:
            print(
                f"WARNING: CSTR solver convergence issue | "
                f"success={result.success}, objective={result.fun:.3e}"
            )

        Fch4_out, Fh2_out = result.x
        Ftot = Fch4_out + Fh2_out
        Pch4 = Fch4_out * press / Ftot
        Ph2 = press - Pch4

        particle_fractions, time_in_reactor = RTD(tres)
        average_a = float(np.sum(deactivation(temp, Pch4, Ph2, time_in_reactor, r_constants, da_constants)*particle_fractions))


        if verbose:
            log(f"Simulation complete: a = {average_a:.4f}, Pch4 = {Pch4 / 1e5:.3f} bar, Ph2 = {Ph2 / 1e5:.3f} bar")

        Fc_out = initial_reaction_rate(temp, Pch4, Ph2, r_constants) * average_a * cat_wt * 1000

        return Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres

    # Mode: constant catalyst rate
    if mode_type == 'constant_rate':
        cat_rate = mode_value
        Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres = simulate(mode_value, Fch4_in, Fh2_in, temp,
                                                                                   press, r_constants, da_constants)

    # Mode: iterate catalyst rate to achieve target average deactivation
    elif mode_type == 'constant_da':
        error = 1
        cat_rate = init_guess
        max_iter = 500 # added finite max number of iterations
        iteration = 0

        while error > 0.001 and iteration < max_iter:
            Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres = simulate(cat_rate, Fch4_in, Fh2_in, temp,
                                                                             press, r_constants, da_constants)


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
           raise RuntimeError("CSTR constant_da did not converge")

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

    def update_deactivation(temp, Pch4, Ph2, a_prev, r_constants, da_constants, dt, initial=False):
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
            a = deactivation(temp, Pch4, Ph2, dt, r_constants, da_constants)
        else:
            k1, k2, k3, k4 = RK4constants(temp, Pch4, Ph2, a_prev, r_constants, da_constants, dt)
            a = a_prev + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        return max(min(a, 1), 0.001)

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

        a = update_deactivation(temp, Pch4, Ph2, a_prev, r_constants, da_constants, dt, initial=(i == 0))

        r = initial_reaction_rate(temp, Pch4, Ph2, r_constants) * a

        # Residuals from molar balances
        R_ch4 = Fch4_in - Fch4_out - r * cat_wt * 1000/ n
        R_h2 = Fh2_in - Fh2_out + 2 * r * cat_wt * 1000 / n

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

        # Pressure at the first step
        Ftot_in_first = Fch4_in + Fh2_in

        Ptot_in_first = press
        Pch4_in_first = Fch4_in * press / Ftot_in_first
        Ph2_in_first = Fh2_in * press / Ftot_in_first


        for i in range(n):
            result = minimize(fun=reactor_equations,
                              x0=init_value,
                              args=(Fch4_in, Fh2_in, temp, press, r_constants, da_constants, cat_wt, n, dt, a_prev, i),
                              bounds=[(1e-11, 100000), (1e-11, 100000)])

            if not result.success or result.fun > 1e-8:
                print(
                    f"WARNING: PFR solver convergence issue at segment {i} | "
                    f"success={result.success}, objective={result.fun:.3e}"
                )

            Fch4_out, Fh2_out = result.x
            Ftot = Fch4_out + Fh2_out
            Pch4 = Fch4_out * press / Ftot
            Ph2 = press - Pch4

            k_diag, Keq_diag = reaction_adsorption_constant(temp, r_constants)
            P0 = 100000
            driving_force = Pch4 - (Ph2 ** 2) / (Keq_diag * P0)
            r0_diag = k_diag * driving_force

            if driving_force < 0 or r0_diag < 0:
                print(
                    f"WARNING: Negative reaction driving force at PFR segment {i} | "
                    f"driving_force={driving_force:.3e} Pa, r0={r0_diag:.3e} kmol/(kg_cat s)"
                )

            a[i] = update_deactivation(temp, Pch4, Ph2, a_prev, r_constants, da_constants, dt, initial=(i == 0))
            a_prev = a[i]
            r = initial_reaction_rate(temp, Pch4, Ph2, r_constants) * a[i]
            Fc_out = Fc_in + r * cat_wt * 1000 / n

            # Update inlet flows and pressure for the next segment
            Fch4_in, Fh2_in, Fc_in = Fch4_out, Fh2_out, Fc_out
            press -= p_drop * 1e5 / n  #from bar to Pascal

            if verbose and i % 100 == 0:
                log(f"[Step {i}/{n}] a = {a[i]:.4f}, Pch4 = {Pch4 / 1e5:.3f} bar, Ph2 = {Ph2 / 1e5:.3f} bar")

        return Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, float(np.average(a)), tres


    # Mode: constant catalyst rate
    if mode_type == 'constant_rate':
        cat_rate = mode_value
        Fch4_out, Fh2_out, Fc_out, Pch4, Ph2, average_a, tres = simulate(mode_value, Fch4_in, Fh2_in, temp,
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
        "Pch4": Pch4 / 1e5,  # Convert from Pascal to bar
        "Ph2": Ph2 / 1e5,  # Convert from Pascal to bar
        "Fch4_out": Fch4_out,
        "Fh2_out": Fh2_out,
        "Fc_out": Fc_out,
        "average_a": average_a,
        "cat_rate": cat_rate,
        "tres": tres
    }
