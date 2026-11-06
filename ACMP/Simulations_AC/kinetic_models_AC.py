import math
import numpy as np
from scipy.optimize import minimize, brentq

def Arrhenius(A, E, T):
    R = 8.3144598e-3  # kJ/(mol K)
    return A * np.exp(-E / (R * T))


def reaction_adsorption_constant(T, r_constants):
    T_k = T + 273.15
    return Arrhenius(*r_constants['k'], T_k)


def equilibrium_constant(T, equilibrium_settings):
    T_k = T + 273.15
    R = 8.3144598  # J/(mol K)

    c0, c1, c2, c3 = equilibrium_settings['delta_g_coefficients']
    delta_g = c0 + c1 * T_k + c2 * T_k**2 + c3 / T_k

    return np.exp(-delta_g / (R * T_k))


def calculate_equilibrium_conversion(T, total_pressure, Fch4_in, Fh2_in, equilibrium_settings):
    if Fch4_in <= 0:
        raise ValueError("The reactor inlet CH4 flow rate must be greater than zero.")

    Kp = equilibrium_constant(T, equilibrium_settings)
    standard_pressure = equilibrium_settings['standard_pressure']

    def equilibrium_residual(conversion):
        Fch4_eq = Fch4_in * (1 - conversion)
        Fh2_eq = Fh2_in + 2 * Fch4_in * conversion
        Ftot_eq = Fch4_eq + Fh2_eq

        Pch4_eq = total_pressure * Fch4_eq / Ftot_eq
        Ph2_eq = total_pressure * Fh2_eq / Ftot_eq
        Q = (Ph2_eq / standard_pressure)**2 / (Pch4_eq / standard_pressure)

        return Q - Kp

    conversion_min = 0.0
    conversion_max = 1.0 - 1e-12
    residual_min = equilibrium_residual(conversion_min)
    residual_max = equilibrium_residual(conversion_max)

    if residual_min >= 0:
        return 0.0

    if residual_max <= 0:
        return conversion_max

    return brentq(equilibrium_residual, conversion_min, conversion_max, xtol=1e-12, rtol=1e-10)


def deactivation_constant(T, da_constants):
    T_k = T + 273.15
    T_c = 190.4
    m = 112 * (T_k / T_c)**(-2.9551)

    Ad_atm, Ed = da_constants['kd']
    atm_to_pa = 101325.0
    Ad_pa = Ad_atm / atm_to_pa**m

    return Arrhenius(Ad_pa, Ed, T_k)


def residenceTime(cat_holdup, cat_rate):
    return cat_holdup / cat_rate


def deactivation(T, Pch4, t, da_constants):
    kd = deactivation_constant(T, da_constants)
    t = np.array(t)

    d = 0.5
    T_k = T + 273.15
    T_c = 190.4
    m = 112 * (T_k / T_c)**(-2.9551)

    denominator = 1 - (1 - d) * kd * Pch4**m * t
    a = np.zeros_like(denominator)

    valid = denominator > 0
    a[valid] = denominator[valid]**(1 / (1 - d))
    a = np.clip(a, 0, 1)

    return a if t.size > 1 else a.item()


def RK4constants(T, Pch4, a_prev, da_constants, dt):
    kd = deactivation_constant(T, da_constants)

    d = 0.5
    T_k = T + 273.15
    T_c = 190.4
    m = 112 * (T_k / T_c)**(-2.9551)

    deactivation_term = kd * Pch4**m

    if deactivation_term <= 0:
        return 0.0, 0.0, 0.0, 0.0

    a_prev = np.clip(a_prev, 0, 1)
    t0 = (1 - a_prev**(1 - d)) / ((1 - d) * deactivation_term)

    def activity_from_time(t):
        denominator = 1 - (1 - d) * deactivation_term * t

        if denominator <= 0:
            return 0.0

        return denominator**(1 / (1 - d))

    def dadt_from_time(t):
        a = activity_from_time(t)
        return -deactivation_term * a**d

    k1 = dt * dadt_from_time(t0)
    k2 = dt * dadt_from_time(t0 + dt / 2)
    k3 = dt * dadt_from_time(t0 + dt / 2)
    k4 = dt * dadt_from_time(t0 + dt)

    return k1, k2, k3, k4


def initial_reaction_rate(T, Pch4, r_constants, equilibrium_factor=1.0):
    k = reaction_adsorption_constant(T, r_constants)
    n = 0.5

    return equilibrium_factor * k * Pch4**n


# CSTR model: separate loops with Brent root finding and warm start
def isothermal_cstr(sim_config, Fch4_in, Fh2_in, p_drop, cat_wt, init_guess, verbose=False, equilibrium_factor_guess=1.0):
    temp = sim_config['reactor']['temperature']
    press = (sim_config['comp1_pressure'] - p_drop) * 1e5
    mode_type = sim_config['mode']['type']
    mode_value = sim_config['mode']['value']
    step_size = sim_config['mode']['step_size']

    r_constants = sim_config['reaction_constants']
    da_constants = sim_config['deactivation_constants']

    Fch4_reactor_in = Fch4_in
    Fh2_reactor_in = Fh2_in
    activity_tolerance = 0.001

    equilibrium_settings = sim_config.get('equilibrium', {})
    equilibrium_enabled = equilibrium_settings.get('enabled', False)
    conversion_tolerance = equilibrium_settings.get('conversion_tolerance', 1e-3)

    def F(t, tau):
        return 1 - math.exp(-t / tau)

    def RTD(tau):
        fractions = []
        times = []
        dt = 0

        while dt < 7 * tau:
            fractions.append(F(dt + tau / 10, tau) - F(dt, tau))
            times.append(dt + tau / 20)
            dt += tau / 10

        fractions = np.array(fractions)
        fractions /= np.sum(fractions)
        return fractions, np.array(times)

    def reactor_equations(flows, tres, equilibrium_factor):
        Fch4_out, Fh2_out = flows
        Ftot = Fch4_out + Fh2_out
        Pch4 = Fch4_out * press / Ftot
        particle_fractions, time_in_reactor = RTD(tres)
        average_a = np.sum(deactivation(temp, Pch4, time_in_reactor, da_constants) * particle_fractions)

        r0 = initial_reaction_rate(temp, Pch4, r_constants, equilibrium_factor)
        r = r0 * average_a

        R_ch4 = Fch4_reactor_in - Fch4_out - r * cat_wt * 1000
        R_h2 = Fh2_reactor_in - Fh2_out + 2 * r * cat_wt * 1000

        return R_ch4**2 + R_h2**2

    previous_flows = None

    def simulate(cat_rate, equilibrium_factor):
        nonlocal previous_flows
        tres = residenceTime(cat_wt, cat_rate)

        if previous_flows is None:
            initial_flows = np.array([Fch4_reactor_in, Fh2_reactor_in])
        else:
            initial_flows = previous_flows.copy()

        result = minimize(reactor_equations, x0=initial_flows, args=(tres, equilibrium_factor), bounds=[(1e-11, 1e6), (1e-11, 1e6)])

        # Fallback only if the standard solve does not satisfy the balances.
        if not result.success or result.fun > 1e-8:
            print(f"WARNING: CSTR solver convergence issue | " f"success={result.success}, objective={result.fun:.3e}. " f"Trying alternative initial guesses.")

            fallback_guesses = []

            for X_guess in (0.50, 0.80, 0.95):
                fallback_guesses.append(np.array([Fch4_reactor_in * (1 - X_guess), Fh2_reactor_in + 2 * Fch4_reactor_in * X_guess]))

            if equilibrium_enabled and X_eq is not None:
                X_guess = min(max(X_eq, 0.0), 0.999999)

                fallback_guesses.append(np.array([Fch4_reactor_in * (1 - X_guess), Fh2_reactor_in + 2 * Fch4_reactor_in * X_guess]))

            best_result = result

            for fallback_guess in fallback_guesses:
                fallback_result = minimize(reactor_equations, x0=fallback_guess, args=(tres, equilibrium_factor), bounds=[(1e-11, 1e6), (1e-11, 1e6)])

                if fallback_result.fun < best_result.fun:
                    best_result = fallback_result

                if fallback_result.success and fallback_result.fun <= 1e-8:
                    best_result = fallback_result
                    break

            result = best_result

            if not result.success or result.fun > 1e-8:
                raise RuntimeError(f"CSTR flow solver failed after fallback guesses: " f"success={result.success}, objective={result.fun:.3e}, " f"CatRate={cat_rate:.6f} kg/s, feq={equilibrium_factor:.6f}")

        Fch4_out, Fh2_out = result.x
        previous_flows = result.x.copy()

        Ftot = Fch4_out + Fh2_out
        Pch4 = Fch4_out * press / Ftot
        Ph2 = press - Pch4

        particle_fractions, time_in_reactor = RTD(tres)
        average_a = float(np.sum(deactivation(temp, Pch4, time_in_reactor, da_constants) * particle_fractions))

        r0 = initial_reaction_rate(temp, Pch4, r_constants, equilibrium_factor)
        Fc_out = r0 * average_a * cat_wt * 1000
        conversion = (Fch4_reactor_in - Fch4_out) / Fch4_reactor_in

        return {
            "Fch4_out": Fch4_out,
            "Fh2_out": Fh2_out,
            "Fc_out": Fc_out,
            "Pch4": Pch4,
            "Ph2": Ph2,
            "average_a": average_a,
            "cat_rate": cat_rate,
            "tres": tres,
            "conversion": conversion
        }

    X_eq = None

    if equilibrium_enabled:
        total_pressure = press

        X_eq = calculate_equilibrium_conversion(T=temp, total_pressure=total_pressure, Fch4_in=Fch4_reactor_in, Fh2_in=Fh2_reactor_in, equilibrium_settings=equilibrium_settings)

    equilibrium_factor = float(np.clip(equilibrium_factor_guess, 1e-8, 1.0))

    if mode_type == 'constant_rate':
        cat_rate = mode_value
        max_joint_iterations = 150

        for iteration in range(max_joint_iterations):
            results = simulate(cat_rate, equilibrium_factor)

            if not equilibrium_enabled:
                conversion_converged = True
            elif results['conversion'] > X_eq + conversion_tolerance:
                conversion_converged = False
            elif abs(results['conversion'] - X_eq) <= conversion_tolerance:
                conversion_converged = True
            elif equilibrium_factor >= 1.0 - 1e-10:
                conversion_converged = True
            else:
                conversion_converged = False

            if conversion_converged:
                break

            damping = 0.5
            ratio = X_eq / max(results['conversion'], 1e-12)
            equilibrium_factor = float(np.clip(equilibrium_factor * ratio**damping, 1e-8, 1.0))

        else:
            raise RuntimeError(f"The CSTR equilibrium iteration did not converge after {max_joint_iterations} iterations.")

    elif mode_type == 'constant_da':
        equilibrium_factor = float(np.clip(equilibrium_factor_guess, 1e-8, 1.0))
        max_equilibrium_iterations = 350
        equilibrium_damping = 0.5

        for equilibrium_iteration in range(max_equilibrium_iterations):
            previous_flows = None

            def activity_residual(cat_rate):
                nonlocal previous_flows
                previous_flows = None

                results_local = simulate(cat_rate, equilibrium_factor)
                residual = results_local['average_a'] - mode_value

                print(f"CSTR equilibrium iteration {equilibrium_iteration}: "
                      f"DA={results_local['average_a']:.5f} (target={mode_value:.5f}) | "
                      f"CatRate={cat_rate:.5f} kg/s")

                return residual

            cat_rate_lower = max(init_guess * 0.1, 1e-8)
            cat_rate_upper = init_guess * 2.0
            residual_lower = activity_residual(cat_rate_lower)
            residual_upper = activity_residual(cat_rate_upper)
            bracket_iteration = 0
            max_bracket_iterations = 50

            while residual_lower * residual_upper > 0 and bracket_iteration < max_bracket_iterations:
                if residual_lower > 0 and residual_upper > 0:
                    cat_rate_lower *= 0.5
                    residual_lower = activity_residual(cat_rate_lower)

                elif residual_lower < 0 and residual_upper < 0:
                    cat_rate_upper *= 2.0
                    residual_upper = activity_residual(cat_rate_upper)

                bracket_iteration += 1

            if residual_lower * residual_upper > 0:
                raise RuntimeError("Unable to bracket the CSTR catalyst rate for the requested activity.")

            cat_rate = brentq(activity_residual, cat_rate_lower, cat_rate_upper, xtol=1e-6, rtol=1e-8, maxiter=100)

            previous_flows = None
            results = simulate(cat_rate, equilibrium_factor)
            activity_error = results['average_a'] - mode_value

            print(f"CSTR equilibrium iteration {equilibrium_iteration} FINAL ACTIVITY: "
                  f"DA={results['average_a']:.5f} (target={mode_value:.5f}) | "
                  f"CatRate={cat_rate:.5f} kg/s")

            if abs(activity_error) > activity_tolerance:
                raise RuntimeError(f"Brentq did not find a valid CSTR activity solution: " f"DA={results['average_a']:.6f}, target={mode_value:.6f}, " f"CatRate={cat_rate:.6f} kg/s, feq={equilibrium_factor:.6f}")

            if not equilibrium_enabled:
                break

            if results['conversion'] > X_eq + conversion_tolerance:
                conversion_converged = False
            elif abs(results['conversion'] - X_eq) <= conversion_tolerance:
                conversion_converged = True
            elif equilibrium_factor >= 1.0 - 1e-10:
                conversion_converged = True
            else:
                conversion_converged = False

            if conversion_converged:
                break

            ratio = X_eq / max(results['conversion'], 1e-12)
            equilibrium_factor = float(np.clip(equilibrium_factor * ratio**equilibrium_damping, 1e-8, 1.0))

        else:
            raise RuntimeError(f"The CSTR equilibrium iteration did not converge after {max_equilibrium_iterations} iterations.")

    else:
        raise ValueError(f"Unsupported mode type: {mode_type}")

    Fch4_out = results['Fch4_out']
    Fh2_out = results['Fh2_out']
    Fc_out = results['Fc_out']
    Pch4 = results['Pch4']
    Ph2 = results['Ph2']
    average_a = results['average_a']
    cat_rate = results['cat_rate']
    tres = results['tres']
    X_python = results['conversion']

    print("\nFINAL PYTHON CSTR")
    print(f"Fch4_in               = {Fch4_reactor_in:.6e} mol/s")
    print(f"Fh2_in                = {Fh2_reactor_in:.6e} mol/s")
    print(f"Fch4_out              = {Fch4_out:.6e} mol/s")
    print(f"Fh2_out               = {Fh2_out:.6e} mol/s")
    print(f"Average activity      = {average_a:.6f}")
    print(f"Catalyst rate         = {cat_rate:.6f} kg/s")
    print(f"Residence time        = {tres:.3f} s")

    return {
        "Pch4": Pch4 / 1e5,
        "Ph2": Ph2 / 1e5,
        "Fch4_out": Fch4_out,
        "Fh2_out": Fh2_out,
        "Fc_out": Fc_out,
        "average_a": average_a,
        "cat_rate": cat_rate,
        "tres": tres,
        "equilibrium_factor": equilibrium_factor,
        "equilibrium_conversion": X_eq
    }


# PFR model: local equilibrium cap applied segment by segment
# The original Abbas kinetic expression is kept unchanged, a local equilibrium factor is introduced only when the kinetic solution exceeds the local equilibrium conversion.
def isothermal_pfr(sim_config, Fch4_in, Fh2_in, p_drop, cat_wt, init_guess, n=501, verbose=False, equilibrium_factor_guess=1.0):
    temp = sim_config['reactor']['temperature']
    press_in = sim_config['comp1_pressure'] * 1e5

    mode_type = sim_config['mode']['type']
    mode_value = sim_config['mode']['value']
    step_size = sim_config['mode']['step_size']

    r_constants = sim_config['reaction_constants']
    da_constants = sim_config['deactivation_constants']

    Fch4_reactor_in = Fch4_in
    Fh2_reactor_in = Fh2_in
    activity_tolerance = 0.001

    equilibrium_settings = sim_config.get('equilibrium', {})
    equilibrium_enabled = equilibrium_settings.get('enabled', False)

    local_equilibrium_tolerance = equilibrium_settings.get('local_equilibrium_tolerance', 1e-6)

    def update_deactivation(Pch4, a_prev, dt, initial=False):
        if initial:
            a = deactivation(temp, Pch4, dt, da_constants)

        else:
            k1, k2, k3, k4 = RK4constants(temp, Pch4, a_prev, da_constants, dt)
            a = a_prev + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        return float(np.clip(a, 0, 1))

    def reactor_equations(flows, Fch4_seg_in, Fh2_seg_in, local_press, dt, a_prev, i):
        Fch4_seg_out, Fh2_seg_out = flows

        Ftot = Fch4_seg_out + Fh2_seg_out
        Pch4 = Fch4_seg_out * local_press / Ftot

        a_local = update_deactivation(Pch4, a_prev, dt, initial=(i == 0))

        r0 = initial_reaction_rate(temp, Pch4, r_constants)
        r = r0 * a_local

        R_ch4 = Fch4_seg_in - Fch4_seg_out - r * cat_wt * 1000 / n
        R_h2 = Fh2_seg_in - Fh2_seg_out + 2 * r * cat_wt * 1000 / n

        return R_ch4**2 + R_h2**2

    def simulate(cat_rate):
        tres = residenceTime(cat_wt, cat_rate)
        dt = tres / n
        a = np.zeros(n)

        Fch4_seg_in = Fch4_reactor_in
        Fh2_seg_in = Fh2_reactor_in
        Fc_in = 0.0

        a_prev = 1.0
        local_press = press_in

        Fch4_out = None
        Fh2_out = None
        Fc_out = None
        Pch4 = None
        Ph2 = None

        capped_segments = 0
        first_capped_segment = None

        local_equilibrium_factors = np.ones(n)
        max_balance_residual = 0.0

        for i in range(n):
            X_eq_local = None

            if equilibrium_enabled:
                X_eq_local = calculate_equilibrium_conversion(T=temp, total_pressure=local_press, Fch4_in=Fch4_seg_in, Fh2_in=Fh2_seg_in, equilibrium_settings=equilibrium_settings)

            init_value = np.array([Fch4_seg_in, Fh2_seg_in])

            max_ch4 = max(Fch4_seg_in * 1.05, 1.0)

            max_h2 = max(Fh2_seg_in + 2 * Fch4_seg_in * 1.05, 1.0)

            result = minimize(reactor_equations, x0=init_value, args=(Fch4_seg_in, Fh2_seg_in, local_press, dt, a_prev, i), bounds=[(1e-11, max_ch4), (1e-11, max_h2)])

            if not result.success or result.fun > 1e-8:
                print(f"WARNING: PFR solver convergence issue at segment {i} | " f"success={result.success}, objective={result.fun:.3e}")

            Fch4_out, Fh2_out = result.x

            X_kin_local = (Fch4_seg_in - Fch4_out) / max(Fch4_seg_in, 1e-12)

            # Equilibrium becomes active only when the kinetic conversion exceeds the local Xeq by more than the dedicated local tolerance.
            equilibrium_limited = (equilibrium_enabled and X_kin_local > X_eq_local + local_equilibrium_tolerance)

            if equilibrium_limited:
                reacted_ch4 = Fch4_seg_in * X_eq_local

                Fch4_out = Fch4_seg_in - reacted_ch4
                Fh2_out = Fh2_seg_in + 2 * reacted_ch4

                capped_segments += 1

                if first_capped_segment is None:
                    first_capped_segment = i

            Ftot = Fch4_out + Fh2_out
            Pch4 = Fch4_out * local_press / Ftot
            Ph2 = local_press - Pch4

            a[i] = update_deactivation(Pch4, a_prev, dt, initial=(i == 0))

            r0_actual = initial_reaction_rate(temp, Pch4, r_constants)

            r_kin_actual = r0_actual * a[i]

            reacted_ch4 = max(Fch4_seg_in - Fch4_out, 0.0)

            segment_cat_term = cat_wt * 1000 / n

            if equilibrium_limited:
                if segment_cat_term <= 0 or r_kin_actual <= 0:
                    raise RuntimeError(f"Invalid equilibrium-limited PFR state at segment {i}: " f"segment_cat_term={segment_cat_term:.6e}, " f"r_kin={r_kin_actual:.6e}.")

                r_required = reacted_ch4 / segment_cat_term

                local_factor_raw = r_required / r_kin_actual

                if local_factor_raw > 1.0 + 1e-6:
                    raise RuntimeError(f"PFR local equilibrium factor exceeds 1 at segment {i}: " f"feq_local={local_factor_raw:.6f}. ")

                local_equilibrium_factors[i] = float(np.clip(local_factor_raw, 1e-12, 1.0))

                r_effective = (local_equilibrium_factors[i] * r_kin_actual)

            else:
                local_equilibrium_factors[i] = 1.0
                r_effective = r_kin_actual

            R_ch4_check = (Fch4_seg_in - Fch4_out - r_effective * segment_cat_term)

            R_h2_check = (Fh2_seg_in - Fh2_out + 2 * r_effective * segment_cat_term)

            segment_balance_residual = max(abs(R_ch4_check), abs(R_h2_check))

            max_balance_residual = max(max_balance_residual, segment_balance_residual)

            if equilibrium_limited:
                balance_tolerance = (1e-8 * max(1.0, abs(Fch4_seg_in), abs(Fh2_seg_in)))

                if segment_balance_residual > balance_tolerance:
                    raise RuntimeError(f"PFR equilibrium-capped balance check failed at segment {i}: " f"R_CH4={R_ch4_check:.6e}, " f"R_H2={R_h2_check:.6e}.")

            a_prev = a[i]

            Fc_out = Fc_in + reacted_ch4

            Fch4_seg_in = Fch4_out
            Fh2_seg_in = Fh2_out
            Fc_in = Fc_out

            local_press -= p_drop * 1e5 / n

            if verbose and i % 100 == 0:
                print(f"[PFR segment {i}/{n}] a={a[i]:.4f} | Pch4={Pch4 / 1e5:.4f} bar | Ph2={Ph2 / 1e5:.4f} bar")

        average_a = float(np.average(a))

        conversion = (Fch4_reactor_in - Fch4_out) / Fch4_reactor_in

        return {
            "Fch4_out": Fch4_out,
            "Fh2_out": Fh2_out,
            "Fc_out": Fc_out,
            "Pch4": Pch4,
            "Ph2": Ph2,
            "average_a": average_a,
            "cat_rate": cat_rate,
            "tres": tres,
            "conversion": conversion,
            "capped_segments": capped_segments,
            "first_capped_segment": first_capped_segment,
            "min_local_equilibrium_factor": float(np.min(local_equilibrium_factors)),
            "max_balance_residual": max_balance_residual
        }

    X_eq = None

    if equilibrium_enabled:
        outlet_pressure = max(press_in - p_drop * 1e5, 1e-12)

        X_eq = calculate_equilibrium_conversion(T=temp, total_pressure=outlet_pressure, Fch4_in=Fch4_reactor_in, Fh2_in=Fh2_reactor_in, equilibrium_settings=equilibrium_settings)

    if mode_type == 'constant_rate':
        cat_rate = mode_value
        results = simulate(cat_rate)

    elif mode_type == 'constant_da':
        cat_rate = init_guess
        max_iterations = 450

        for iteration in range(max_iterations):
            results = simulate(cat_rate)

            activity_error = (results['average_a'] - mode_value)

            print(f"PFR iteration {iteration}: DA={results['average_a']:.5f} "
                  f"(target={mode_value:.5f}) | CatRate={cat_rate:.5f} kg/s")

            if abs(activity_error) <= activity_tolerance:
                break

            if results['average_a'] > mode_value:
                cat_rate -= (abs(activity_error) * init_guess * step_size)

            else:
                cat_rate += (abs(activity_error) * init_guess * step_size)

            if cat_rate <= 0:
                raise RuntimeError("Catalyst rate has reached 0.")

        else:
            raise RuntimeError(f"The PFR activity iteration did not converge after {max_iterations} iterations.")

    else:
        raise ValueError(f"Unsupported mode type: {mode_type}")

    Fch4_out = results['Fch4_out']
    Fh2_out = results['Fh2_out']
    Fc_out = results['Fc_out']
    Pch4 = results['Pch4']
    Ph2 = results['Ph2']
    average_a = results['average_a']
    cat_rate = results['cat_rate']
    tres = results['tres']
    X_python = results['conversion']

    print("\nFINAL PYTHON PFR")
    print(f"Fch4_in               = {Fch4_reactor_in:.6e} mol/s")
    print(f"Fh2_in                = {Fh2_reactor_in:.6e} mol/s")
    print(f"Fch4_out              = {Fch4_out:.6e} mol/s")
    print(f"Fh2_out               = {Fh2_out:.6e} mol/s")
    print(f"Average activity      = {average_a:.6f}")
    print(f"Catalyst rate         = {cat_rate:.6f} kg/s")
    print(f"Residence time        = {tres:.3f} s")

    return {
        "Pch4": Pch4 / 1e5,
        "Ph2": Ph2 / 1e5,
        "Fch4_out": Fch4_out,
        "Fh2_out": Fh2_out,
        "Fc_out": Fc_out,
        "average_a": average_a,
        "cat_rate": cat_rate,
        "tres": tres,
        "conversion": X_python,
        "equilibrium_factor": 1.0,
        "equilibrium_conversion": X_eq,
        "capped_segments": results['capped_segments'],
        "first_capped_segment": results['first_capped_segment'],
        "min_local_equilibrium_factor": results['min_local_equilibrium_factor'],
        "max_balance_residual": results['max_balance_residual']
    }
