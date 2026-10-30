from math import pi
from time import sleep
from pathlib import Path

import numpy as np
import pandas as pd
import win32com.client as win32

from scipy.optimize import minimize, brentq

from kinetic_models_AC import isothermal_cstr, isothermal_pfr



def initialize_aspen(file_path, sim_config):
    aspen = win32.gencache.EnsureDispatch('Apwn.Document')
    aspen.InitFromArchive2(file_path)
    aspen.Visible = sim_config['open_aspen']
    aspen.SuppressDialogs = True

    return aspen


def find_node(aspen, node):
    return aspen.Tree.FindNode('\\Data\\' + node).Value


def assign_node(aspen, node, value):
    aspen.Tree.FindNode('\\Data\\' + node).Value = value


def get_temp_values(reactor_temp, reactor_type):

    cstr_ranges = [
        (800, 25, 80, 635),
        (825, 45, 80, 635),
        (850, 60, 75, 665),
        (875, 65, 70, 705),
        (900, 65, 70, 790),
        (925, 65, 70, 825),
        (950, 65, 70, 895)
    ]

    pfr_ranges = [
        (800, 25, 75, 635),
        (825, 65, 75, 640),
        (850, 65, 75, 725),
        (875, 65, 75, 790),
        (900, 65, 75, 840),
        (925, 65, 75, 875),
        (950, 65, 75, 905)
    ]

    ranges = pfr_ranges if reactor_type == 'PFR' else cstr_ranges

    for max_temp, hx1_value, hx2_temp, hx3_temp in ranges:
        if reactor_temp <= max_temp:
            return hx1_value, hx2_temp, hx3_temp

    raise ValueError("Invalid reactor temperature")


def assign_temperatures(aspen, reactor_temp, heating_method, reactor_type):
    hx1_value, hx2_temp, hx3_temp = get_temp_values(reactor_temp, reactor_type)

    assign_node(aspen, r'Blocks\HX-1\Input\VALUE', hx1_value)

    if heating_method != 'electric':
        assign_node(aspen, r'Blocks\HX-2A\Input\TEMP', hx2_temp)
        assign_node(aspen, r'Blocks\HX-3A\Input\TEMP', hx3_temp)


def kinetics_settings(aspen, sim_config):
    A, E = sim_config['reaction_constants']['k']

    assign_node(aspen, r'Reactions\Reactions\CMP\Input\PRE_EXP\1', A)
    assign_node(aspen, r'Reactions\Reactions\CMP\Input\ACT_ENERGY\1', E)


def reactor_size(sim_config, solid_c_wt):
    def objective(vars, cat_wt, cat_density, solid_c_wt, solid_c_density):
        dia, height = vars

        aspect_ratio = height / dia

        particle_volume = cat_wt / cat_density + solid_c_wt / solid_c_density
        min_vol = particle_volume * 3

        actual_vol = pi * (dia / 2)**2 * height
        cross_area = pi * (dia / 2)**2

        p_drop = ((cat_wt + solid_c_wt) * 9.8067 * 1e-5 * (4 / 3) / cross_area)

        punishment1 = 1000 if actual_vol < min_vol else 0
        punishment2 = 1000 if aspect_ratio < 2.5 else 0

        return abs(sim_config['comp1_pressure'] - p_drop - sim_config['reactor']['pressure']) + punishment1 + punishment2

    num_units = 1

    while True:
        adj_cat_wt = sim_config['catalyst']['weight'] / num_units
        cat_density = sim_config['catalyst']['density']
        adj_solid_c_wt = solid_c_wt / num_units
        solid_c_density = sim_config['solid_c_density']

        result = minimize(objective, x0=np.array([10, 30]), args=(adj_cat_wt, cat_density, adj_solid_c_wt, solid_c_density), bounds=[(0.1, 10), (1, 30)])

        dia, height = result.x

        aspect_ratio = height / dia
        volume = pi * (dia / 2)**2 * height
        cross_area = pi * (dia / 2)**2

        p_drop = ((adj_cat_wt + adj_solid_c_wt) * 9.8067 * 1e-5 * (4 / 3) / cross_area)

        min_volume = (adj_cat_wt / cat_density + adj_solid_c_wt / solid_c_density) * (3 / 1)

        if (min_volume < volume <= 2400 and min_volume <= 2400 and aspect_ratio >= 2.5 and sim_config['comp1_pressure'] - p_drop >= 1.05):
            break

        num_units += 1

    return num_units, volume, p_drop


def initialize_simulation(sim_config):
    path = Path().resolve()

    file_paths = {
        'PFR': {
            'electric': fr'{path}\aspen_plus\AC\pfr-electric-AC.apw',
            'ch4': fr'{path}\aspen_plus\AC\pfr-ch4-AC.apw',
            'h2': fr'{path}\aspen_plus\AC\pfr-h2-AC.apw'
        },
        'CSTR': {
            'electric': fr'{path}\aspen_plus\AC\cstr-electric-AC.apw',
            'ch4': fr'{path}\aspen_plus\AC\cstr-ch4-AC.apw',
            'h2': fr'{path}\aspen_plus\AC\cstr-h2-AC.apw'
        }
    }

    reactor_type = sim_config['reactor']['type']
    reactor_temp = sim_config['reactor']['temperature']
    heating_method = sim_config['reactor']['heating']

    aspen = initialize_aspen(file_paths[reactor_type][heating_method], sim_config)

    assign_node(aspen, r'Setup\Sim-Options\Input\Unit Set', 'PBT1')
    assign_node(aspen, r'Convergence\Conv-Options\Input\TOL', sim_config['tolerance'])

    assign_node(aspen, r'Flowsheeting Options\Design-Spec\CAPACITY\Input\EXPR2', sim_config['plant_capacity'])

    capacity_lower = (sim_config['plant_capacity'] / 100) * (20000 if heating_method == 'h2' else 15000)

    capacity_upper = (sim_config['plant_capacity'] / 100) * (30000 if heating_method == 'h2' else 20000)

    assign_node(aspen, r'Flowsheeting Options\Design-Spec\CAPACITY\Input\LOWER', capacity_lower)
    assign_node(aspen, r'Flowsheeting Options\Design-Spec\CAPACITY\Input\UPPER', capacity_upper)

    assign_node(aspen, r'Blocks\HEATER-1\Input\TEMP', reactor_temp)
    assign_node(aspen, r'Blocks\HEATER-2\Input\TEMP', reactor_temp)

    assign_node(aspen, r'Blocks\COMP-1\Input\PRES', sim_config['comp1_pressure'])
    assign_node(aspen, r'Blocks\COMP-2\Input\PRES', sim_config['comp2_pressure'])

    if heating_method == 'ch4':
        assign_node(aspen, r'Blocks\COMP-4\Input\PRES', sim_config['comp4_pressure'])

    assign_node(aspen, r'Blocks\CYCLONE\Input\DES_PDROP', 0.05)

    if sim_config['mode']['type'] == 'constant_rate':
        cat_rate = sim_config['mode']['value']
    else:
        cat_rate = sim_config['initial_guess']

    assign_node(aspen, r'Streams\CAT-FEED\Input\TOTFLOW\CIPSD', cat_rate)
    assign_node(aspen, r'Streams\CAT-FEED\Input\PRES\CIPSD', sim_config['comp1_pressure'])
    assign_node(aspen, r'Blocks\REACTOR\Input\CATWT', sim_config['catalyst']['weight'])

    assign_temperatures(aspen, reactor_temp, heating_method, reactor_type)

    num_reacs, reac_vol_unit, p_drop = reactor_size(sim_config, 10000)

    print(f"Number of reactors: {num_reacs}, " f"Volume per reactor = {reac_vol_unit:.4g}, " f"Reactor pressure = {(sim_config['comp1_pressure'] - p_drop):.4g}")

    if reactor_type == 'PFR':
        assign_node(aspen, r'Blocks\REACTOR\Input\REAC_TEMP', reactor_temp)
        assign_node(aspen, r'Blocks\REACTOR\Input\PDROP', p_drop)

        reac_dia = ((4 * reac_vol_unit * num_reacs) / (3 * pi))**(1 / 3)
        reac_height = 3 * reac_dia

        assign_node(aspen, r'Blocks\REACTOR\Input\DIAM', reac_dia)
        assign_node(aspen, r'Blocks\REACTOR\Input\LENGTH', reac_height)

    elif reactor_type == 'CSTR':
        assign_node(aspen, r'Blocks\REACTOR\Input\TEMP', reactor_temp)
        assign_node(aspen, r'Blocks\REACTOR\Input\PRES', sim_config['comp1_pressure'] - p_drop)
        assign_node(aspen, r'Blocks\REACTOR\Input\VOL', reac_vol_unit * num_reacs)

    kinetics_settings(aspen, sim_config)

    if heating_method != 'electric':
        if heating_method == 'ch4':
            rankine_flow = 9500 if reactor_type == 'PFR' else 10000
        else:
            rankine_flow = 12500 if reactor_type == 'PFR' else 15000

        assign_node(aspen, r'Blocks\RANKINE\Data\Streams\R2\Input\TOTFLOW\MIXED', sim_config['plant_capacity'] / 100 * rankine_flow)

        assign_node(aspen, r'Blocks\CMBST\Data\Blocks\FURNACE\Input\TEMP', reactor_temp + 20)

    return aspen, num_reacs, p_drop


def run_aspen_simulation(aspen, max_time=100):
    aspen.Reinit()
    aspen.Engine.Run2(1)

    elapsed_time = 0

    while aspen.Engine.IsRunning == 1 and elapsed_time < max_time:
        sleep(0.5)
        elapsed_time += 0.5

    if aspen.Engine.IsRunning == 1:
        aspen.Engine.Stop()


def get_flow_rates(aspen):
    Fch4_in = find_node(aspen, r'Streams\PS6\Output\MOLEFLOW\MIXED\CH4') * 1000 / 3600
    Fh2_in = find_node(aspen, r'Streams\PS6\Output\MOLEFLOW\MIXED\H2') * 1000 / 3600

    Mc_out = find_node(aspen, r'Streams\SOLID-C\Output\MASSFLOW\CIPSD\C')

    Fch4_out = find_node(aspen, r'Streams\PS7\Output\MOLEFLOW\MIXED\CH4') * 1000 / 3600
    Fh2_out = find_node(aspen, r'Streams\PS7\Output\MOLEFLOW\MIXED\H2') * 1000 / 3600

    return Fch4_in, Fh2_in, Mc_out, Fch4_out, Fh2_out


def update_aspen(sim_config, aspen, num_reacs, reac_vol_unit, p_drop, cat_rate_unit, deactivation, equilibrium_factor=1.0):

    A, _ = sim_config['reaction_constants']['k']

    assign_node(aspen, r'Reactions\Reactions\CMP\Input\PRE_EXP\1', A * deactivation * equilibrium_factor)

    assign_node(aspen, r'Streams\CAT-FEED\Input\TOTFLOW\CIPSD', cat_rate_unit * num_reacs)

    if sim_config['reactor']['type'] == 'PFR':
        assign_node(aspen, r'Blocks\REACTOR\Input\PDROP', p_drop)

        reac_dia = ((4 * reac_vol_unit * num_reacs) / (3 * pi))**(1 / 3)
        reac_height = 3 * reac_dia

        assign_node(aspen, r'Blocks\REACTOR\Input\DIAM', reac_dia)
        assign_node(aspen, r'Blocks\REACTOR\Input\LENGTH', reac_height)

    else:
        assign_node(aspen, r'Blocks\REACTOR\Input\PRES', sim_config['comp1_pressure'] - p_drop)
        assign_node(aspen, r'Blocks\REACTOR\Input\VOL', reac_vol_unit * num_reacs)


def run_reactor(aspen, sim_config, num_reacs, p_drop):
    Fch4_in, Fh2_in, Mc, Fch4_out, Fh2_out = get_flow_rates(aspen)

    Fch4_unit = Fch4_in / num_reacs
    Fh2_unit = Fh2_in / num_reacs

    cat_wt_unit = sim_config['catalyst']['weight'] / num_reacs
    cat_rate_unit = sim_config['initial_guess'] / num_reacs

    reactor_type = sim_config['reactor']['type']
    equilibrium_settings = sim_config.get('equilibrium', {})
    equilibrium_enabled = equilibrium_settings.get('enabled', False)

    aspen_conversion_tolerance = equilibrium_settings.get(
        'aspen_conversion_tolerance',
        equilibrium_settings.get('conversion_tolerance', 1e-3)
    )
    aspen_matching_factor_min = equilibrium_settings.get(
        'aspen_matching_factor_min',
        1e-8
    )
    aspen_matching_factor_max = equilibrium_settings.get(
        'aspen_matching_factor_max',
        5.0
    )
    aspen_matching_growth = equilibrium_settings.get(
        'aspen_matching_growth',
        1.5
    )

    equilibrium_factor_guess = 1.0
    aspen_matching_factor = 1.0

    if reactor_type == 'PFR':
        results = isothermal_pfr(
            sim_config,
            Fch4_unit,
            Fh2_unit,
            p_drop,
            cat_wt_unit,
            cat_rate_unit,
            equilibrium_factor_guess=1.0
        )
    else:
        results = isothermal_cstr(
            sim_config,
            Fch4_unit,
            Fh2_unit,
            p_drop,
            cat_wt_unit,
            cat_rate_unit,
            equilibrium_factor_guess=equilibrium_factor_guess
        )

    deactivation = results['average_a']
    cat_rate = results['cat_rate'] * num_reacs
    tres = results['tres']

    if reactor_type == 'CSTR':
        equilibrium_factor_guess = results.get('equilibrium_factor', 1.0)

    error = 10.0
    aspen_match_error = np.inf if reactor_type == 'PFR' and equilibrium_enabled else 0.0
    iteration = 0
    max_iterations = 20

    Fch4_out_unit = Fch4_out / num_reacs
    Fh2_out_unit = Fh2_out / num_reacs


    Fch4_total = Fch4_in
    Fh2_total = Fh2_in

    while iteration < max_iterations:
        solid_c_wt = (Mc / 3600) * tres

        num_reacs, reac_vol_unit, p_drop = reactor_size(sim_config, solid_c_wt)

        cat_wt_unit = sim_config['catalyst']['weight'] / num_reacs
        cat_rate_unit = cat_rate / num_reacs

        if reactor_type == 'PFR':
            X_python_target = results['conversion']

            def evaluate_aspen_factor(factor):
                update_aspen(
                    sim_config,
                    aspen,
                    num_reacs,
                    reac_vol_unit,
                    p_drop,
                    cat_rate_unit,
                    deactivation,
                    factor
                )

                run_aspen_simulation(aspen)

                flow_values = get_flow_rates(aspen)
                Fch4_eval, _, _, Fch4_out_eval, _ = flow_values
                X_aspen_eval = (
                    Fch4_eval - Fch4_out_eval
                ) / max(Fch4_eval, 1e-12)

                return X_aspen_eval - X_python_target, X_aspen_eval, flow_values

            if not equilibrium_enabled:
                aspen_matching_factor = 1.0
                _, X_aspen, flow_values = evaluate_aspen_factor(aspen_matching_factor)
                aspen_match_error = 0.0

            else:
                factor_guess = float(np.clip(
                    aspen_matching_factor,
                    aspen_matching_factor_min,
                    aspen_matching_factor_max
                ))

                residual_guess, X_aspen_guess, flow_values_guess = evaluate_aspen_factor(factor_guess)

                best_factor = factor_guess
                best_residual = residual_guess
                best_X_aspen = X_aspen_guess
                best_flow_values = flow_values_guess

                if abs(residual_guess) <= aspen_conversion_tolerance:
                    aspen_matching_factor = factor_guess
                    X_aspen = X_aspen_guess
                    flow_values = flow_values_guess

                elif residual_guess > 0:
                    # Aspen conversion is above the Python target: search below the current factor.
                    factor_high = factor_guess
                    factor_low = factor_guess
                    bracket_found = False

                    while factor_low > aspen_matching_factor_min * (1 + 1e-12):
                        factor_low = max(
                            factor_low / aspen_matching_growth,
                            aspen_matching_factor_min
                        )

                        residual_low, X_aspen_low, flow_values_low = evaluate_aspen_factor(factor_low)

                        if abs(residual_low) < abs(best_residual):
                            best_factor = factor_low
                            best_residual = residual_low
                            best_X_aspen = X_aspen_low
                            best_flow_values = flow_values_low

                        if residual_low <= 0:
                            bracket_found = True
                            break

                        factor_high = factor_low

                    if bracket_found:
                        def conversion_residual(factor):
                            residual, _, _ = evaluate_aspen_factor(factor)
                            return residual

                        aspen_matching_factor = brentq(
                            conversion_residual,
                            factor_low,
                            factor_guess,
                            xtol=1e-6,
                            rtol=1e-6,
                            maxiter=30
                        )

                        _, X_aspen, flow_values = evaluate_aspen_factor(aspen_matching_factor)
                    else:
                        aspen_matching_factor = best_factor
                        X_aspen = best_X_aspen
                        flow_values = best_flow_values
                        print(
                            "Aspen PFR conversion remains above the Python target "
                            "at the minimum allowed Aspen matching factor."
                        )

                else:

                    factor_low = factor_guess
                    factor_high = factor_guess
                    bracket_found = False

                    while factor_high < aspen_matching_factor_max * (1 - 1e-12):
                        factor_high = min(
                            factor_high * aspen_matching_growth,
                            aspen_matching_factor_max
                        )

                        residual_high, X_aspen_high, flow_values_high = evaluate_aspen_factor(factor_high)

                        if abs(residual_high) < abs(best_residual):
                            best_factor = factor_high
                            best_residual = residual_high
                            best_X_aspen = X_aspen_high
                            best_flow_values = flow_values_high

                        if residual_high >= 0:
                            bracket_found = True
                            break

                    if bracket_found:
                        def conversion_residual(factor):
                            residual, _, _ = evaluate_aspen_factor(factor)
                            return residual

                        aspen_matching_factor = brentq(
                            conversion_residual,
                            factor_guess,
                            factor_high,
                            xtol=1e-6,
                            rtol=1e-6,
                            maxiter=30
                        )

                        _, X_aspen, flow_values = evaluate_aspen_factor(aspen_matching_factor)
                    else:
                        aspen_matching_factor = best_factor
                        X_aspen = best_X_aspen
                        flow_values = best_flow_values
                        print(
                            "Aspen PFR conversion remains below the Python target "
                            "at the maximum allowed Aspen matching factor."
                        )

                aspen_match_error = abs(X_aspen - X_python_target)

            Fch4, Fh2, Mc, Fch4_out, Fh2_out = flow_values

            print(
                f"Aspen PFR match: X_python={X_python_target:.6f} | "
                f"X_aspen={X_aspen:.6f} | "
                f"f_aspen={aspen_matching_factor:.6f} | "
                f"match error={aspen_match_error * 100:.6f} percentage points"
            )

        else:
            update_aspen(
                sim_config,
                aspen,
                num_reacs,
                reac_vol_unit,
                p_drop,
                cat_rate_unit,
                deactivation,
                equilibrium_factor_guess
            )

            run_aspen_simulation(aspen)
            Fch4, Fh2, Mc, Fch4_out, Fh2_out = get_flow_rates(aspen)

        Fch4_unit_new = Fch4 / num_reacs
        Fh2_unit_new = Fh2 / num_reacs

        Fch4_out_unit = Fch4_out / num_reacs
        Fh2_out_unit = Fh2_out / num_reacs

        error = max(
            abs((Fch4_total - Fch4) / Fch4_total),
            abs((Fh2_total - Fh2) / Fh2_total)
        )

        if sim_config['mode']['type'] == 'constant_rate':
            if reactor_type == 'PFR':
                print(
                    f"Iteration {iteration}: Error = {error * 100:.4g}%, "
                    f"Deactivation Factor = {deactivation:.4g}, "
                    f"Aspen factor = {aspen_matching_factor:.6f}"
                )
            else:
                print(
                    f"Iteration {iteration}: Error = {error * 100:.4g}%, "
                    f"Deactivation Factor = {deactivation:.4g}, "
                    f"Equilibrium factor = {equilibrium_factor_guess:.6f}"
                )
        else:
            if reactor_type == 'PFR':
                print(
                    f"Iteration {iteration}: Error = {error * 100:.4g}%, "
                    f"Catalyst rate = {cat_rate:.4g}, "
                    f"Aspen factor = {aspen_matching_factor:.6f}"
                )
            else:
                print(
                    f"Iteration {iteration}: Error = {error * 100:.4g}%, "
                    f"Catalyst rate = {cat_rate:.4g}, "
                    f"Equilibrium factor = {equilibrium_factor_guess:.6f}"
                )

        print(
            f"Number of reactors: {num_reacs}, "
            f"Volume per reactor = {reac_vol_unit:.4g}, "
            f"Reactor pressure = {(sim_config['comp1_pressure'] - p_drop):.4g}"
        )

        Fch4_total = Fch4
        Fh2_total = Fh2
        Fch4_unit = Fch4_unit_new
        Fh2_unit = Fh2_unit_new

        flow_converged = error < 0.001
        match_converged = (
            reactor_type != 'PFR'
            or not equilibrium_enabled
            or aspen_match_error <= aspen_conversion_tolerance
        )

        if flow_converged and match_converged:
            break

        cat_rate_unit = cat_rate / num_reacs

        if reactor_type == 'PFR':
            results = isothermal_pfr(
                sim_config,
                Fch4_unit,
                Fh2_unit,
                p_drop,
                cat_wt_unit,
                cat_rate_unit,
                equilibrium_factor_guess=1.0
            )
        else:
            results = isothermal_cstr(
                sim_config,
                Fch4_unit,
                Fh2_unit,
                p_drop,
                cat_wt_unit,
                cat_rate_unit,
                equilibrium_factor_guess=equilibrium_factor_guess
            )

        deactivation = results['average_a']
        cat_rate = results['cat_rate'] * num_reacs
        tres = results['tres']

        if reactor_type == 'CSTR':
            equilibrium_factor_guess = results.get('equilibrium_factor', 1.0)

        iteration += 1

    if iteration >= max_iterations:
        print(
            "WARNING: Aspen-Python outer loop reached the maximum number "
            "of iterations before full convergence."
        )

    cat_rate_unit = cat_rate / num_reacs

    if reactor_type == 'CSTR':
        results = isothermal_cstr(
            sim_config,
            Fch4_unit,
            Fh2_unit,
            p_drop,
            cat_wt_unit,
            cat_rate_unit,
            equilibrium_factor_guess=equilibrium_factor_guess
        )

    # Use the last Python solution consistently in the returned values.
    deactivation = results['average_a']
    cat_rate = results['cat_rate'] * num_reacs
    tres = results['tres']

    print()
    print('CH4 reactor outlet rel. error: ', (Fch4_out_unit - results['Fch4_out']) * 100 / results['Fch4_out'], '%')
    print('H2 reactor outlet rel. error: ', (Fh2_out_unit - results['Fh2_out']) * 100 / results['Fh2_out'], '%')

    return deactivation, cat_rate, tres, num_reacs, reac_vol_unit, p_drop

def collect_results(sim_config, aspen, deactivation, tres, num_reacs, reac_vol, arrays):
    run_status = find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR')

    error_message = ''

    if run_status == 0:
        arrays['CH4 Feed [kg/h]'][0] = find_node(aspen, r'Streams\CH4-FEED\Output\MASSFLMX_GAS')
        arrays['H2 Prod [kg/h]'][0] = find_node(aspen, r'Streams\H2-PROD\Output\MASSFLMX_GAS')
        arrays['C Prod [kg/h]'][0] = find_node(aspen, r'Streams\SOLID-C\Output\MASSFLOW\CIPSD\C')
        arrays['CH4 Recycle [kg/h]'][0] = find_node(aspen, r'Streams\PS15\Output\MASSFLMX_GAS')

        arrays['CH4 Reactor Feed [kmol/h]'][0] = find_node(aspen, r'Streams\PS6\Output\MOLEFLOW\MIXED\CH4')

        arrays['CH4 Reactor Outlet [kmol/h]'][0] = find_node(aspen, r'Streams\PS7\Output\MOLEFLOW\MIXED\CH4')

        arrays['Cat Rate [kg/h]'][0] = find_node(aspen, r'Streams\CAT-FEED\Output\MASSFLMX\CIPSD')

        arrays['Cat Res. Time [s]'][0] = tres
        arrays['Reactor Deactivation'][0] = deactivation
        arrays['Number of Reactors'][0] = num_reacs
        arrays['Reactor Volume [m3]'][0] = reac_vol

        arrays['Comp-1 Power [kW]'][0] = find_node(aspen, r'Blocks\COMP-1\Output\WNET')
        arrays['Comp-2 Power [kW]'][0] = find_node(aspen, r'Blocks\COMP-2\Output\WNET')

        arrays['HeatX-1 Area [m2]'][0] = find_node(aspen, r'Blocks\HX-1\Output\HX_AREAP')

        arrays['Cooler-1 Area [m2]'][0] = find_node(aspen, r'Blocks\COOLER-1\Output\HX_AREAP')
        arrays['Cooler-2 Area [m2]'][0] = find_node(aspen, r'Blocks\COOLER-2\Output\HX_AREAP')

        arrays['Cool Water 1 [kJ/s]'][0] = find_node(aspen, r'Blocks\COOLER-1\Output\HX_DUTY')
        arrays['Cool Water 2 [kJ/s]'][0] = find_node(aspen, r'Blocks\COOLER-2\Output\HX_DUTY')

        arrays['Heater-1 Duty [MW]'][0] = find_node(aspen, r'Blocks\HEATER-1\Output\QNET') / 1000
        arrays['Heater-2 Duty [MW]'][0] = find_node(aspen, r'Blocks\HEATER-2\Output\QNET') / 1000
        arrays['Reactor Duty [MW]'][0] = find_node(aspen, r'Blocks\REACTOR\Output\QCALC') / 1000

        arrays['Cyclone Volumetric [m3/s]'][0] = find_node(aspen, r'Streams\PS7\Output\VOLFLMX2')

        arrays['PSA Mole Feed [kmol/hr]'][0] = find_node(aspen, r'Blocks\PSA\Output\BAL_MOLI_TFL')

        arrays['Error Status'][0] = run_status

        if sim_config['reactor']['heating'] != 'electric':
            arrays['Comp-3 Power [kW]'][0] = find_node(aspen, r'Blocks\CMBST\Data\Blocks\COMP-3\Output\WNET')

            arrays['Comp-3 Vol. Flow [m3/h]'][0] = (find_node(aspen, r'Blocks\CMBST\Data\Blocks\COMP-3\Output\VFLOW') * 3600)

            arrays['Fuel Feed [kg/h]'][0] = find_node(aspen, r'Streams\FUELFEED\Output\MASSFLMX_GAS')

            arrays['CO2 Emission [kg/h]'][0] = find_node(aspen, r'Streams\FLUEGAS\Output\MASSFLOW\MIXED\CO2')

            arrays['Cool Water 3 [kJ/s]'][0] = find_node(aspen, r'Blocks\RANKINE\Data\Blocks\COOLER-3\Output\HX_DUTY')

            arrays['HeatX-2 Area [m2]'][0] = find_node(aspen, r'Blocks\HEATX-2\Data\Blocks\HX-2\Output\HX_AREAP')

            arrays['HeatX-3 Area [m2]'][0] = find_node(aspen, r'Blocks\HEATX-3\Data\Blocks\HX-3\Output\HX_AREAP')

            arrays['Cooler-3 Area [m2]'][0] = find_node(aspen, r'Blocks\RANKINE\Data\Blocks\COOLER-3\Output\HX_AREAP')

            arrays['Turb-1 Power [kW]'][0] = find_node(aspen, r'Blocks\RANKINE\Data\Blocks\TURB-1\Output\WNET')

            arrays['Pump Flow [L/s]'][0] = (find_node(aspen, r'Blocks\RANKINE\Data\Blocks\PUMP\Output\VFLOW') * 1000)

            cw_tot = (find_node(aspen, r'Blocks\COOLER-1\Output\HX_DUTY') + find_node(aspen, r'Blocks\COOLER-2\Output\HX_DUTY') - find_node(aspen, r'Blocks\COMP-2\Output\DUTY_OUT') + find_node(aspen, r'Blocks\RANKINE\Data\Blocks\COOLER-3\Output\HX_DUTY'))

            if sim_config['reactor']['heating'] == 'ch4':
                arrays['Comp-4 Power [kW]'][0] = find_node(aspen, r'Blocks\COMP-4\Output\WNET')
                arrays['Cool Water 4 [kJ/s]'][0] = find_node(aspen, r'Blocks\COOLER-4\Output\HX_DUTY')
                arrays['Cooler-4 Area [m2]'][0] = find_node(aspen, r'Blocks\COOLER-4\Output\HX_AREAP')

                arrays['CO2 with CCS [kg/h]'][0] = find_node(aspen, r'Streams\CLEANGAS\Output\MASSFLOW\MIXED\CO2')

                cw_tot += (find_node(aspen, r'Blocks\COOLER-4\Output\HX_DUTY') - find_node(aspen, r'Blocks\COMP-4\Output\DUTY_OUT'))

            else:
                arrays['Comp-4 Power [kW]'][0] = 0
                arrays['Cool Water 4 [kJ/s]'][0] = 0
                arrays['Cooler-4 Area [m2]'][0] = 0
                arrays['CO2 with CCS [kg/h]'][0] = 0

            arrays['Cool Water Total [kJ/s]'][0] = cw_tot

        else:
            arrays['Comp-3 Power [kW]'][0] = 0
            arrays['Comp-3 Vol. Flow [m3/h]'][0] = 0
            arrays['Comp-4 Power [kW]'][0] = 0

            arrays['Fuel Feed [kg/h]'][0] = 0
            arrays['CO2 Emission [kg/h]'][0] = 0
            arrays['CO2 with CCS [kg/h]'][0] = 0

            arrays['Cool Water 3 [kJ/s]'][0] = 0
            arrays['Cool Water 4 [kJ/s]'][0] = 0

            arrays['HeatX-2 Area [m2]'][0] = 0
            arrays['HeatX-3 Area [m2]'][0] = 0
            arrays['Cooler-3 Area [m2]'][0] = 0
            arrays['Cooler-4 Area [m2]'][0] = 0

            arrays['Turb-1 Power [kW]'][0] = 0
            arrays['Pump Flow [L/s]'][0] = 0

            arrays['Cool Water Total [kJ/s]'][0] = (find_node(aspen, r'Blocks\COOLER-1\Output\HX_DUTY') + find_node(aspen, r'Blocks\COOLER-2\Output\HX_DUTY') - find_node(aspen, r'Blocks\COMP-2\Output\DUTY_OUT'))

    else:
        for key in arrays:
            if key == 'Error Status':
                arrays[key][0] = run_status

                error_message = str(find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR\1') + " " + find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR\2') + " " + find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR\3') + " " + find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR\4'))

            else:
                arrays[key][0] = float('inf')

    return run_status, error_message


def run_scenario(sim_config):
    aspen, num_reacs, p_drop = initialize_simulation(sim_config)

    run_aspen_simulation(aspen)

    deactivation, cat_rate, tres, num_reacs, reac_vol_unit, p_drop = run_reactor(aspen, sim_config, num_reacs, p_drop)

    result_keys = [
        'CH4 Feed [kg/h]',
        'H2 Prod [kg/h]',
        'C Prod [kg/h]',
        'CH4 Recycle [kg/h]',
        'CH4 Reactor Feed [kmol/h]',
        'CH4 Reactor Outlet [kmol/h]',
        'Cat Rate [kg/h]',
        'Cat Res. Time [s]',
        'Reactor Deactivation',
        'Number of Reactors',
        'Reactor Volume [m3]',
        'Comp-1 Power [kW]',
        'Comp-2 Power [kW]',
        'Comp-3 Power [kW]',
        'Comp-3 Vol. Flow [m3/h]',
        'Comp-4 Power [kW]',
        'HeatX-1 Area [m2]',
        'HeatX-2 Area [m2]',
        'HeatX-3 Area [m2]',
        'Turb-1 Power [kW]',
        'Pump Flow [L/s]',
        'Cool Water Total [kJ/s]',
        'Cooler-1 Area [m2]',
        'Cooler-2 Area [m2]',
        'Cooler-3 Area [m2]',
        'Cooler-4 Area [m2]',
        'Cool Water 1 [kJ/s]',
        'Cool Water 2 [kJ/s]',
        'Cool Water 3 [kJ/s]',
        'Cool Water 4 [kJ/s]',
        'Fuel Feed [kg/h]',
        'CO2 Emission [kg/h]',
        'CO2 with CCS [kg/h]',
        'Heater-1 Duty [MW]',
        'Heater-2 Duty [MW]',
        'Reactor Duty [MW]',
        'Cyclone Volumetric [m3/s]',
        'PSA Mole Feed [kmol/hr]',
        'Error Status'
    ]

    arrays = {key: np.empty(1) for key in result_keys}

    run_status, error_message = collect_results(sim_config, aspen, deactivation, tres, num_reacs, reac_vol_unit, arrays)

    if error_message:
        print(f"Error message: {error_message}")

    sleep(1.5)

    if run_status == 0:
        aspen.Reinit()
        aspen.Application.Quit()

    reactor_type = sim_config['reactor']['type']
    capacity = sim_config['plant_capacity']
    heat_source = sim_config['reactor']['heating']
    temp = sim_config['reactor']['temperature']
    press = sim_config['reactor']['pressure']
    catwt = sim_config['catalyst']['weight']
    psa_press = sim_config['comp2_pressure']

    if sim_config['mode']['type'] == 'constant_rate':
        file_name = (f'sim-{reactor_type}-{capacity}TPD-{heat_source}-' f't{temp}-p{round(press, 3)}-c{round(catwt, 1)}-' f'cr{round(cat_rate, 3)}-psa{psa_press}.csv')

    else:
        file_name = (f'sim-{reactor_type}-{capacity}TPD-{heat_source}-' f't{temp}-p{round(press, 3)}-c{round(catwt, 1)}-' f'da{round(deactivation, 2)}-psa{psa_press}.csv')

    df = pd.DataFrame({key: arrays[key] for key in result_keys})
    df.to_csv(file_name, mode='w', index=False)

    return deactivation, cat_rate, tres
