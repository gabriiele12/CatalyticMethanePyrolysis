from math import pi, log
from time import sleep
import win32com.client as win32
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize
import numpy as np
from kinetic_models_Ni import isothermal_cstr, isothermal_pfr


def initialize_aspen(file_path, sim_config):
    """
    Initialize an Aspen simulation from a specified file path and configuration.

    Parameters:
    file_path (str): The path to the Aspen simulation file.
    sim_config (dict): A dictionary containing simulation configuration options.

    Returns:
    aspen: An initialized Aspen simulation object.
    """
    aspen = win32.gencache.EnsureDispatch('Apwn.Document')
    aspen.InitFromArchive2(file_path)
    aspen.Visible = sim_config['open_aspen']
    aspen.SuppressDialogs = True  # Suppress windows dialogs
    return aspen


def find_node(aspen, node):
    """
    Find and return the value of a specified node in the Aspen simulation tree.

    Parameters:
    aspen: The Aspen simulation object.
    node (str): The path to the node within the Aspen simulation tree.

    Returns:
    The value of the specified node.
    """
    return aspen.Tree.FindNode('\\Data\\' + node).Value


def assign_node(aspen, node, value):
    """
    Assign a value to a specified node in the Aspen simulation tree.

    Parameters:
    aspen: The Aspen simulation object.
    node (str): The path to the node within the Aspen simulation tree.
    value: The value to assign to the specified node.
    """
    aspen.Tree.FindNode('\\Data\\' + node).Value = value


def initialize_simulation(sim_config):
    """
    Initialize an Aspen simulation based on the provided simulation configuration.

    Parameters:
    sim_config (dict): A dictionary containing simulation configuration options.

    Returns:
    tuple: A tuple containing the initialized Aspen simulation object,
           the number of reactors, and the pressure drop.
    """
    path = Path().resolve()
    file_paths = {
        'PFR': {
            'electric': fr'{path}\aspen_plus\pfr-electric-Ni.apw',
            'ch4': fr'{path}\aspen_plus\pfr-ch4-Ni.apw',
            'h2': fr'{path}\aspen_plus\pfr-h2-Ni.apw'
        },
        'CSTR': {
            'electric': fr'{path}\aspen_plus\cstr-electric-Ni.apw',
            'ch4': fr'{path}\aspen_plus\cstr-ch4-Ni.apw',
            'h2': fr'{path}\aspen_plus\cstr-h2-Ni.apw'
        }
    }

    reactor_type = sim_config['reactor']['type']
    reactor_temp = sim_config['reactor']['temperature']
    heating_method = sim_config['reactor']['heating']
    aspen = initialize_aspen(file_paths[reactor_type][heating_method], sim_config)

    # Set general simulation options
    assign_node(aspen, r'Setup\Sim-Options\Input\Unit Set', 'PBT1')
    assign_node(aspen, r'Convergence\Conv-Options\Input\TOL', sim_config['tolerance'])
    assign_node(aspen, r'Flowsheeting Options\Design-Spec\CAPACITY\Input\EXPR2', sim_config['plant_capacity'])
    capacity_lower = (sim_config['plant_capacity'] / 100) * (15000 if heating_method == 'h2' else 16000)
    capacity_upper = (sim_config['plant_capacity'] / 100) * (30000 if heating_method == 'h2' else 20000)
    assign_node(aspen, r'Flowsheeting Options\Design-Spec\CAPACITY\Input\LOWER', capacity_lower)
    assign_node(aspen, r'Flowsheeting Options\Design-Spec\CAPACITY\Input\UPPER', capacity_upper)


    # Set reactor and heater temperatures
    for node in [r'Blocks\HEATER-1\Input\TEMP', r'Blocks\HEATER-2\Input\TEMP']:
        assign_node(aspen, node, reactor_temp)

    # Set pressure values
    if heating_method == 'ch4':
        assign_node(aspen, fr'Blocks\COMP-1\Input\PRES', sim_config[f'comp1_pressure'])
        assign_node(aspen, fr'Blocks\COMP-2\Input\PRES', sim_config[f'comp2_pressure'])
        assign_node(aspen, fr'Blocks\COMP-4\Input\PRES', sim_config[f'comp4_pressure'])
    else:
        assign_node(aspen, fr'Blocks\COMP-1\Input\PRES', sim_config[f'comp1_pressure'])
        assign_node(aspen, fr'Blocks\COMP-2\Input\PRES', sim_config[f'comp2_pressure'])
    assign_node(aspen, r'Blocks\CYCLONE\Input\DES_PDROP', 0.05)

    # Set catalyst feed properties
    mode_type = sim_config['mode']['type']
    if mode_type == 'constant_rate':
        cat_rate = sim_config['mode']['value']
    elif mode_type == 'constant_da':
        cat_rate = sim_config['initial_guess']
    else:
        raise ValueError(f"Unsupported mode type: {mode_type}")

    assign_node(aspen, r'Streams\CAT-FEED\Input\TOTFLOW\CIPSD', cat_rate) #assigned in kg/s
    assign_node(aspen, r'Streams\CAT-FEED\Input\PRES\CIPSD', sim_config['comp1_pressure'])
    assign_node(aspen, r'Blocks\REACTOR\Input\CATWT', sim_config['catalyst']['weight'])

    assign_temperatures(aspen, reactor_temp, heating_method, reactor_type)

    # Reactor sizing
    num_reacs, reac_vol_unit, p_drop = reactor_size(sim_config, 10000)
    print(f"Number of reactors: {num_reacs}, Volume per reactor = {reac_vol_unit:.4g}, "
          f"Reactor pressure = {(sim_config['comp1_pressure'] - p_drop):.4g}")

    # Set reactor-specific settings
    if reactor_type == 'PFR':
        assign_node(aspen, r'Blocks\REACTOR\Input\REAC_TEMP', reactor_temp)
        assign_node(aspen, r'Blocks\REACTOR\Input\PDROP', p_drop)
        reac_dia = ((4 * reac_vol_unit * num_reacs) / (3 * pi)) ** (1 / 3)
        reac_height = 3 * reac_dia
        assign_node(aspen, r'Blocks\REACTOR\Input\DIAM', reac_dia)
        assign_node(aspen, r'Blocks\REACTOR\Input\LENGTH', reac_height)
    elif reactor_type == 'CSTR':
        assign_node(aspen, r'Blocks\REACTOR\Input\TEMP', reactor_temp)
        assign_node(aspen, r'Blocks\REACTOR\Input\PRES', sim_config['comp1_pressure'] - p_drop)
        assign_node(aspen, r'Blocks\REACTOR\Input\VOL', reac_vol_unit * num_reacs)
    else:
        raise ValueError(f"Unsupported reactor type: {reactor_type}")

    kinetics_settings(aspen, sim_config)

    if heating_method != 'electric':
        if heating_method == 'ch4':
            if reactor_type == 'PFR':
                assign_node(aspen, r'Blocks\RANKINE\Data\Streams\R2\Input\TOTFLOW\MIXED',
                            (sim_config['plant_capacity'] / 100) * 9500)
            else:
                assign_node(aspen, r'Blocks\RANKINE\Data\Streams\R2\Input\TOTFLOW\MIXED',
                              (sim_config['plant_capacity'] / 100) * 10000)
        elif heating_method == 'h2':
            if reactor_type == 'PFR':
                assign_node(aspen, r'Blocks\RANKINE\Data\Streams\R2\Input\TOTFLOW\MIXED',
                            (sim_config['plant_capacity'] / 100) * 12500)
            else:
                assign_node(aspen, r'Blocks\RANKINE\Data\Streams\R2\Input\TOTFLOW\MIXED',
                              (sim_config['plant_capacity'] / 100) * 15000)
        assign_node(aspen, r'Blocks\CMBST\Data\Blocks\FURNACE\Input\TEMP', reactor_temp+20)

    return aspen, num_reacs, p_drop

def get_temp_values(reactor_temp, reactor_type):
    cstr_ranges = [
        (575, 45, 100, 220),
        (600, 45, 100, 200),
        (625, 45, 95, 190),
        (650, 50, 95, 200),
        (675, 50, 90, 220),
        (700, 50, 90, 240)
    ]

    pfr_ranges = [
        (575, 45, 100, 220),
        (600, 50, 100, 180),
        (625, 45, 90, 190),
        (650, 50, 85, 230),
        (675, 50, 80, 245),
        (700, 60, 80, 280)
    ]

    if reactor_type == 'PFR':
        for max_temp, hx1_value, hx2_temp, hx3_temp in pfr_ranges:
            if reactor_temp <= max_temp:
                return hx1_value, hx2_temp, hx3_temp
    elif reactor_type == 'CSTR':
        for max_temp, hx1_value, hx2_temp, hx3_temp in cstr_ranges:
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
    """
    Set the reaction constants in the Aspen simulation based on the provided simulation configuration.

    Parameters:
    aspen: The Aspen simulation object.
    sim_config (dict): A dictionary containing simulation configuration options, including reaction constants.
    """
    reaction_constants = sim_config['reaction_constants']

    # Set reaction constants for k
    A, E = reaction_constants['k']
    assign_node(aspen, r'Reactions\Reactions\MP\Input\PRE_EXP\1', A/1000)
    assign_node(aspen, r'Reactions\Reactions\MP\Input\ACT_ENERGY\1', E)

    # Set reaction constants for K_p
    Ap, Ep = reaction_constants['K_p']
    assign_node(aspen, r'Reactions\Reactions\MP\Input\DF2_A\1', log(1 / Ap))
    assign_node(aspen, r'Reactions\Reactions\MP\Input\DF2_B\1', Ep * 1000 / 8.3144598)

    # Set reaction constants for K_ch4
    A_ch4, E_ch4 = reaction_constants['K_ch4']
    assign_node(aspen, r'Reactions\Reactions\MP\Input\ADS_A\1\3', log(A_ch4))
    assign_node(aspen, r'Reactions\Reactions\MP\Input\ADS_B\1\3', -E_ch4 * 1000 / 8.3144598)

    # Set reaction constants for K_h2
    A_h2, E_h2 = reaction_constants['K_h2']
    assign_node(aspen, r'Reactions\Reactions\MP\Input\ADS_A\1\2', log(A_h2))
    assign_node(aspen, r'Reactions\Reactions\MP\Input\ADS_B\1\2', -E_h2 * 1000 / 8.3144598)


def reactor_size(sim_config, solid_c_wt):
    """
    Calculate the size and number of reactors required based on the simulation configuration and solid catalyst weight.

    Parameters:
    sim_config (dict): A dictionary containing simulation configuration options.
    solid_c_wt (float): The weight of the solid catalyst.

    Returns:
    tuple: A tuple containing the number of reactors, the volume of each reactor, and the pressure drop.
    """
    def objective(vars, cat_wt, cat_density, solid_c_wt, solid_c_density):
        """
        Objective function to minimize for reactor sizing.

        Parameters:
        vars (tuple): A tuple containing the diameter and height of the reactor.
        cat_wt (float): The weight of the catalyst.
        cat_density (float): The density of the catalyst.
        solid_c_wt (float): The weight of the solid catalyst.
        solid_c_density (float): The density of the solid catalyst.

        Returns:
        float: The calculated objective value.
        """
        dia, height = vars
        aspect_ratio = height / dia
        particle_volume = (cat_wt / cat_density) + (solid_c_wt / solid_c_density)
        min_vol = particle_volume * 3
        actual_vol = pi * (dia / 2) ** 2 * height
        cross_area = pi * (dia / 2) ** 2
        p_drop = (cat_wt + solid_c_wt) * 9.8067 * 1e-5 * (4/3) / cross_area
        punishment1 = 1000 if actual_vol < min_vol else 0
        punishment2 = 1000 if aspect_ratio < 2.5 else 0
        return (abs((sim_config['comp1_pressure'] - p_drop - sim_config['reactor']['pressure']))
                + punishment1 + punishment2)

    num_units = 1

    while True:
        adj_cat_wt = sim_config['catalyst']['weight'] / num_units
        cat_density = sim_config['catalyst']['density']
        adj_solid_c_wt = solid_c_wt / num_units
        solid_c_density = sim_config['solid_c_density']

        result = minimize(fun=objective,
                           x0=np.array([10, 30]),
                          args=(adj_cat_wt, cat_density, adj_solid_c_wt, solid_c_density),
                          bounds=[(0.1, 10), (1, 30)])

        dia, height = result.x[0], result.x[1]
        aspect_ratio = height / dia
        volume = pi * (dia / 2) ** 2 * height
        cross_area = pi * (dia / 2) ** 2
        p_drop = (adj_cat_wt + adj_solid_c_wt) * 9.8067 * 1e-5 * (4/3) / cross_area
        min_volume = ((adj_cat_wt / cat_density) + (adj_solid_c_wt / solid_c_density)) * (3 / 1)

        if (min_volume < volume <= 2400 and
                min_volume <= 2400 and
                aspect_ratio >= 2.5 and
                sim_config['comp1_pressure'] - p_drop >= 1.05):
            break

        num_units += 1

    return num_units, volume, p_drop


def run_aspen_simulation(aspen, max_time=100):
    """
    Run the Aspen simulation for a specified maximum time.

    Parameters:
    aspen: The Aspen simulation object.
    max_time (int): The maximum time to run the simulation in seconds. Default is 100 seconds.
    """
    aspen.Reinit()
    aspen.Engine.Run2(1)
    elapsed_time = 0

    # Loop until the simulation is complete or the maximum time is reached
    while aspen.Engine.IsRunning == 1 and elapsed_time < max_time:
        sleep(0.5)
        elapsed_time += 0.5

    # Stop the simulation if it is still running after the maximum time
    if aspen.Engine.IsRunning == 1:
        aspen.Engine.Stop()


def get_flow_rates(aspen):
    """
    Retrieve the flow rates of various components from the Aspen simulation.

    Parameters:
    aspen: The Aspen simulation object.

    Returns:
    tuple: A tuple containing the flow rates of CH4 and H2 at the reactor feed and product streams,
           and the mass flow rate of solid carbon.
    """
    Fch4_in = find_node(aspen, r'Streams\PS6\Output\MOLEFLOW\MIXED\CH4') * 1000 / 3600
    Fh2_in = find_node(aspen, r'Streams\PS6\Output\MOLEFLOW\MIXED\H2') * 1000 / 3600
    Mc_out = find_node(aspen, r'Streams\SOLID-C\Output\MASSFLOW\CIPSD\C')
    Fch4_out = find_node(aspen, r'Streams\PS7\Output\MOLEFLOW\MIXED\CH4') * 1000 / 3600
    Fh2_out = find_node(aspen, r'Streams\PS7\Output\MOLEFLOW\MIXED\H2') * 1000 / 3600
    return Fch4_in, Fh2_in, Mc_out, Fch4_out, Fh2_out


def update_aspen(sim_config, aspen, num_reacs, reac_vol_unit, p_drop, cat_rate_unit, deactivation):
    """
    Update the Aspen simulation with the provided configuration and parameters.

    Parameters:
    sim_config (dict): A dictionary containing simulation configuration options.
    aspen: The Aspen simulation object.
    num_reacs (int): The number of reactors.
    reac_vol_unit (float): The volume of each reactor unit.
    p_drop (float): The pressure drop across the reactor.
    cat_rate_unit (float): The catalyst feed rate per unit.
    deactivation (float): The deactivation factor for the catalyst.

    Updates the Aspen simulation with the specified parameters for either PFR or CSTR reactor types.
    """
    A, _ = sim_config['reaction_constants']['k']
    assign_node(aspen, r'Reactions\Reactions\MP\Input\PRE_EXP\1', A * deactivation / 1000)
    assign_node(aspen, r'Streams\CAT-FEED\Input\TOTFLOW\CIPSD', cat_rate_unit * num_reacs) #assigned in kg/s

    if sim_config['reactor']['type'] == 'PFR':
        # PFR-specific settings
        assign_node(aspen, r'Blocks\REACTOR\Input\PDROP', p_drop)
        reac_dia = ((4 * reac_vol_unit * num_reacs) / (3 * pi)) ** (1 / 3)
        reac_height = 3 * reac_dia
        assign_node(aspen, r'Blocks\REACTOR\Input\DIAM', reac_dia)
        assign_node(aspen, r'Blocks\REACTOR\Input\LENGTH', reac_height)

    elif sim_config['reactor']['type'] == 'CSTR':
        # CSTR-specific settings
        assign_node(aspen, r'Blocks\REACTOR\Input\PRES', sim_config['comp1_pressure'] - p_drop)
        assign_node(aspen, r'Blocks\REACTOR\Input\VOL', reac_vol_unit * num_reacs)


def run_reactor(aspen, sim_config, num_reacs, p_drop):
    """
    Run the reactor simulation and iteratively update the Aspen simulation until convergence.

    Parameters:
    aspen: The Aspen simulation object.
    sim_config (dict): A dictionary containing simulation configuration options.
    num_reacs (int): The number of reactors.
    p_drop (float): The pressure drop across the reactor.

    Returns:
    tuple: A tuple containing the deactivation factor, catalyst rate, residence time,
           number of reactors, reactor volume per unit, and pressure drop.
    """
    Fch4_in, Fh2_in, Mc, _, _ = get_flow_rates(aspen)

    cat_wt_unit = sim_config['catalyst']['weight'] / num_reacs
    Fch4_unit = Fch4_in / num_reacs
    Fh2_unit = Fh2_in / num_reacs
    cat_rate_unit = sim_config['initial_guess'] / num_reacs

    reactor_func = isothermal_pfr if sim_config['reactor']['type'] == 'PFR' else isothermal_cstr
    results = reactor_func(sim_config, Fch4_unit, Fh2_unit, p_drop, cat_wt_unit, cat_rate_unit)
    deactivation, cat_rate, tres = results['average_a'], results['cat_rate'] * num_reacs, results['tres']

    # Initialize error parameters for the iterative process
    error, max_iterations, iteration = 10, 20, 0
    Fch4_total = Fch4_in
    Fh2_total = Fh2_in

    # Iteratively run simulations until convergence or maximum iterations reached
    while error > 0.001 and iteration < max_iterations:
        solid_c_wt = (Mc / 3600) * tres

        num_reacs, reac_vol_unit, p_drop = reactor_size(sim_config, solid_c_wt)
        cat_rate_unit = cat_rate / num_reacs

        update_aspen(sim_config, aspen, num_reacs, reac_vol_unit, p_drop, cat_rate_unit, deactivation)

        run_aspen_simulation(aspen)

        Fch4, Fh2, Mc, Fch4_out, Fh2_out = get_flow_rates(aspen)

        Fch4_unit_new = Fch4 / num_reacs
        Fh2_unit_new = Fh2 / num_reacs
        Fch4_out_unit = Fch4_out / num_reacs
        Fh2_out_unit = Fh2_out / num_reacs
        cat_wt_unit = sim_config['catalyst']['weight'] / num_reacs

        # Calculate the error between old and new total flow rates
        error = max(abs((Fch4_total - Fch4) / Fch4_total),
                    abs((Fh2_total - Fh2) / Fh2_total))
        if sim_config['mode']['type'] == 'constant_rate':
            print(fr"Iteration {iteration}: Error = {error*100:.4g}%, Deactivation Factor = {deactivation:.4g}")
        elif sim_config['mode']['type'] == 'constant_da':
            print(fr"Iteration {iteration}: Error = {error*100:.4g}%, Catalyst rate = {cat_rate:.4g}")
        print(f"Number of reactors: {num_reacs}, Volume per reactor = {reac_vol_unit:.4g}, "
              f"Reactor pressure = {(sim_config['comp1_pressure'] - p_drop):.4g}")

        Fch4_total, Fh2_total = Fch4, Fh2
        Fch4_unit, Fh2_unit = Fch4_unit_new, Fh2_unit_new

        if error < 0.001:
            break

        results = reactor_func(sim_config, Fch4_unit, Fh2_unit, p_drop, cat_wt_unit, cat_rate_unit)
        deactivation, cat_rate, tres = results['average_a'], results['cat_rate'] * num_reacs, results['tres']

        iteration += 1

    results = reactor_func(sim_config, Fch4_unit, Fh2_unit, p_drop, cat_wt_unit, cat_rate_unit)

    print('\n')
    # print(Fch4_out_unit, Fh2_out_unit)
    # print(results['Fch4_out'], results['Fh2_out'])
    print('CH4 reactor outlet rel. error: ', (Fch4_out_unit - results['Fch4_out'])*100/results['Fch4_out'],'%')
    print('H2 reactor outlet rel. error: ', (Fh2_out_unit - results['Fh2_out'])*100/results['Fh2_out'],'%')

    return deactivation, cat_rate, tres, num_reacs, reac_vol_unit, p_drop


def collect_results(sim_config, aspen, deactivation, tres, num_reacs, reac_vol, arrays):
    """
    Collect results from the Aspen simulation and store them in the provided arrays.

    Parameters:
    sim_config (dict): A dictionary containing simulation configuration options.
    aspen: The Aspen simulation object.
    deactivation (float): The deactivation factor for the catalyst.
    tres (float): The residence time of the catalyst.
    num_reacs (int): The number of reactors.
    reac_vol (float): The volume of each reactor.
    arrays (dict): A dictionary to store the collected results.

    Returns:
    tuple: A tuple containing the run status and an error message (if any).
    """
    run_status = find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR')
    error_message = ''
    if run_status == 0:
        arrays['CH4 Feed [kg/h]'][0] = find_node(aspen, r'Streams\CH4-FEED\Output\MASSFLMX_GAS')  # kg/h
        arrays['H2 Prod [kg/h]'][0] = find_node(aspen, r'Streams\H2-PROD\Output\MASSFLMX_GAS')  # kg/h
        arrays['C Prod [kg/h]'][0] = find_node(aspen, r'Streams\SOLID-C\Output\MASSFLOW\CIPSD\C')  # kg/h: MASSFLMX\CIPSD  MASSFLOW\CIPSD\C
        arrays['CH4 Recycle [kg/h]'][0] = find_node(aspen, r'Streams\PS15\Output\MASSFLMX_GAS')  # kg/h
        arrays['CH4 Reactor Feed [kmol/h]'][0] = find_node(aspen, r'Streams\PS6\Output\MOLEFLOW\MIXED\CH4')  # kmol/h
        arrays['CH4 Reactor Outlet [kmol/h]'][0] = find_node(aspen, r'Streams\PS7\Output\MOLEFLOW\MIXED\CH4')  # kmol/h
        arrays['Cat Rate [kg/h]'][0] = find_node(aspen, r'Streams\CAT-FEED\Output\MASSFLMX\CIPSD')  # read in kg/h

        arrays['Cat Res. Time [s]'][0] = tres
        arrays['Reactor Deactivation'][0] = deactivation
        arrays['Number of Reactors'][0] = num_reacs
        arrays['Reactor Volume [m3]'][0] = reac_vol

        arrays['Comp-1 Power [kW]'][0] = find_node(aspen, r'Blocks\COMP-1\Output\WNET')  # kW
        arrays['Comp-2 Power [kW]'][0] = find_node(aspen, r'Blocks\COMP-2\Output\WNET')  # kW
        arrays['HeatX-1 Area [m2]'][0] = find_node(aspen, r'Blocks\HX-1\Output\HX_AREAP')  # m2
        arrays['Cooler-1 Area [m2]'][0] = find_node(aspen, r'Blocks\COOLER-1\Output\HX_AREAP')  # m2
        arrays['Cooler-2 Area [m2]'][0] = find_node(aspen, r'Blocks\COOLER-2\Output\HX_AREAP')  # m2
        arrays['Cool Water 1 [kJ/s]'][0] = find_node(aspen, r'Blocks\COOLER-1\Output\HX_DUTY')  # kJ/s
        arrays['Cool Water 2 [kJ/s]'][0] = find_node(aspen, r'Blocks\COOLER-2\Output\HX_DUTY')  # kJ/s
        arrays['Heater-1 Duty [MW]'][0] = find_node(aspen, r'Blocks\HEATER-1\Output\QNET') / 1000  # MW
        arrays['Heater-2 Duty [MW]'][0] = find_node(aspen, r'Blocks\HEATER-2\Output\QNET') / 1000  # MW
        arrays['Reactor Duty [MW]'][0] = find_node(aspen, r'Blocks\REACTOR\Output\QCALC') / 1000  # MW
        arrays['Cyclone Volumetric [m3/s]'][0] = find_node(aspen, r'Streams\PS7\Output\VOLFLMX2')  # m3/s
        arrays['PSA Mole Feed [kmol/hr]'][0] = find_node(aspen, r'Blocks\PSA\Output\BAL_MOLI_TFL')  # kmol/hr
        arrays['Error Status'][0] = run_status

        if sim_config['reactor']['heating'] != 'electric':
            arrays['Comp-3 Power [kW]'][0] = find_node(aspen, r'Blocks\CMBST\Data\Blocks\COMP-3\Output\WNET')  # kW
            arrays['Comp-3 Vol. Flow [m3/h]'][0] = find_node(aspen, r'Blocks\CMBST\Data\Blocks\COMP-3\Output\VFLOW') * 3600  # m3/h
            arrays['Fuel Feed [kg/h]'][0] = find_node(aspen, r'Streams\FUELFEED\Output\MASSFLMX_GAS')  # kg/h
            arrays['CO2 Emission [kg/h]'][0] = find_node(aspen, r'Streams\FLUEGAS\Output\MASSFLOW\MIXED\CO2')  # kg/h
            arrays['Cool Water 3 [kJ/s]'][0] = find_node(aspen, r'Blocks\RANKINE\Data\Blocks\COOLER-3\Output\HX_DUTY')  # kJ/s
            arrays['HeatX-2 Area [m2]'][0] = find_node(aspen, r'Blocks\HEATX-2\Data\Blocks\HX-2\Output\HX_AREAP')  # m2
            arrays['HeatX-3 Area [m2]'][0] = find_node(aspen, r'Blocks\HEATX-3\Data\Blocks\HX-3\Output\HX_AREAP')  # m2
            arrays['Cooler-3 Area [m2]'][0] = find_node(aspen, r'Blocks\RANKINE\Data\Blocks\COOLER-3\Output\HX_AREAP')  # m2
            arrays['Turb-1 Power [kW]'][0] = find_node(aspen, r'Blocks\RANKINE\Data\Blocks\TURB-1\Output\WNET')  # kW
            arrays['Pump Flow [L/s]'][0] = find_node(aspen, r'Blocks\RANKINE\Data\Blocks\PUMP\Output\VFLOW') * 1000  # L/s
            cw_tot = (find_node(aspen, r'Blocks\COOLER-1\Output\HX_DUTY') +
                      find_node(aspen, r'Blocks\COOLER-2\Output\HX_DUTY') +
                      find_node(aspen, r'Blocks\RANKINE\Data\Blocks\COOLER-3\Output\HX_DUTY') +
                      find_node(aspen, r'Blocks\COMP-2\Output\DUTY_OUT') * -1)  # kJ/s
            if sim_config['reactor']['heating'] == 'ch4':
                arrays['Comp-4 Power [kW]'][0] = find_node(aspen, r'Blocks\COMP-4\Output\WNET')  # kW
                arrays['Cool Water 4 [kJ/s]'][0] = find_node(aspen, r'Blocks\COOLER-4\Output\HX_DUTY')  # kJ/s
                arrays['Cooler-4 Area [m2]'][0] = find_node(aspen, r'Blocks\COOLER-4\Output\HX_AREAP')  # m2
                arrays['CO2 with CCS [kg/h]'][0] = find_node(aspen,
                                                             r'Streams\CLEANGAS\Output\MASSFLOW\MIXED\CO2')
                cw_tot += (find_node(aspen, r'Blocks\COOLER-4\Output\HX_DUTY') +
                           find_node(aspen, r'Blocks\COMP-4\Output\DUTY_OUT') * -1)
                arrays['Cool Water Total [kJ/s]'][0] = cw_tot  # kg/h
            else:
                arrays['Comp-4 Power [kW]'][0] = 0
                arrays['Cool Water 4 [kJ/s]'][0] = 0
                arrays['Cooler-4 Area [m2]'][0] = 0
                arrays['CO2 with CCS [kg/h]'][0] = 0
                arrays['Cool Water Total [kJ/s]'][0] = cw_tot  # kg/h
        else:
            arrays['Comp-3 Power [kW]'][0] = 0
            arrays['Comp-3 Vol. Flow [m3/h]'][0] = 0
            arrays['Comp-4 Power [kW]'][0] = 0
            arrays['Fuel Feed [kg/h]'][0] = 0
            arrays['CO2 Emission [kg/h]'][0] = 0
            arrays['CO2 with CCS [kg/h]'][0] = 0
            cw_tot = (find_node(aspen, r'Blocks\COOLER-1\Output\HX_DUTY') +
                      find_node(aspen, r'Blocks\COOLER-2\Output\HX_DUTY') +
                      find_node(aspen, r'Blocks\COMP-2\Output\DUTY_OUT')*-1)  # kJ/s
            arrays['Cool Water 3 [kJ/s]'][0] = 0
            arrays['Cool Water 4 [kJ/s]'][0] = 0
            arrays['Cool Water Total [kJ/s]'][0] = cw_tot
            arrays['HeatX-2 Area [m2]'][0] = 0  # m2
            arrays['HeatX-3 Area [m2]'][0] = 0  # m2

            arrays['Cooler-3 Area [m2]'][0] = 0
            arrays['Cooler-4 Area [m2]'][0] = 0
            arrays['Turb-1 Power [kW]'][0] = 0
            arrays['Pump Flow [L/s]'][0] = 0

    else:
        for key in arrays:
            if key == 'Error Status':
                arrays[key][0] = run_status
                if run_status != 0:
                    error_message = str(
                    find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR\1') + " " +
                    find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR\2') + " " +
                    find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR\3') + " " +
                    find_node(aspen, r'Results Summary\Run-Status\Output\PER_ERROR\4')
                    )

            else:
                arrays[key][0] = float('inf')

    return run_status, error_message


def run_scenario(sim_config):
    """
    Run a simulation scenario based on the provided configuration.

    Parameters:
    sim_config (dict): A dictionary containing simulation configuration options.

    Returns:
    tuple: A tuple containing the deactivation factor, catalyst rate, and residence time.
    """
    deactivation, cat_rate, tres = 0, 0, 0

    # Initialize the Aspen simulation
    aspen, num_reacs, p_drop = initialize_simulation(sim_config)

    # Run the initial Aspen simulation
    run_aspen_simulation(aspen)

    # Run the reactor simulation and update Aspen until convergence
    deactivation, cat_rate, tres, num_reacs, reac_vol_unit, p_drop = run_reactor(aspen, sim_config, num_reacs, p_drop)

    # Define the keys for the result arrays
    result_keys = ['CH4 Feed [kg/h]', 'H2 Prod [kg/h]', 'C Prod [kg/h]',
                   'CH4 Recycle [kg/h]', 'CH4 Reactor Feed [kmol/h]', 'CH4 Reactor Outlet [kmol/h]',
                   'Cat Rate [kg/h]', 'Cat Res. Time [s]', 'Reactor Deactivation',
                   'Number of Reactors', 'Reactor Volume [m3]',
                   'Comp-1 Power [kW]', 'Comp-2 Power [kW]', 'Comp-3 Power [kW]', 'Comp-3 Vol. Flow [m3/h]',
                   'Comp-4 Power [kW]', 'HeatX-1 Area [m2]', 'HeatX-2 Area [m2]', 'HeatX-3 Area [m2]',
                   'Turb-1 Power [kW]', 'Pump Flow [L/s]', 'Cool Water Total [kJ/s]',
                   'Cooler-1 Area [m2]', 'Cooler-2 Area [m2]', 'Cooler-3 Area [m2]', 'Cooler-4 Area [m2]',
                   'Cool Water 1 [kJ/s]', 'Cool Water 2 [kJ/s]', 'Cool Water 3 [kJ/s]', 'Cool Water 4 [kJ/s]',
                   'Fuel Feed [kg/h]', 'CO2 Emission [kg/h]', 'CO2 with CCS [kg/h]',
                   'Heater-1 Duty [MW]', 'Heater-2 Duty [MW]', 'Reactor Duty [MW]',
                   'Cyclone Volumetric [m3/s]', 'PSA Mole Feed [kmol/hr]',
                   'Error Status']

    # Initialize the result arrays
    arrays = {key: np.empty(1) for key in result_keys}

    # Collect results from the Aspen simulation
    run_status, error_message = collect_results(sim_config, aspen, deactivation, tres, num_reacs, reac_vol_unit, arrays)
    if error_message != '':
        print(f"Error message: {error_message}")
    sleep(1.5)

    # Reinitialize and quit Aspen if the run was successful  #comment to keep open
    if run_status == 0:
        aspen.Reinit()
        aspen.Application.Quit()

    # Generate the file name for saving results
    file_name = ''
    type = sim_config['reactor']['type']
    capacity = sim_config['plant_capacity']
    heat_source = sim_config['reactor']['heating']
    temp = sim_config['reactor']['temperature']
    press = sim_config['reactor']['pressure']
    catwt = sim_config['catalyst']['weight']
    psa_press = sim_config['comp2_pressure']
    # multiplier = sim_config['multiplier']
    if sim_config['mode']['type'] == 'constant_rate':
        file_name = (f'sim-{type}-{capacity}TPD-{heat_source}-t{temp}-p{round(press, 3)}'
                     f'-c{round(catwt, 1)}-cr{round(cat_rate, 3)}-psa{psa_press}.csv')
    elif sim_config['mode']['type'] == 'constant_da':
        file_name = (f'sim-{type}-{capacity}TPD-{heat_source}-t{temp}-p{round(press, 3)}'
                     f'-c{round(catwt, 1)}-da{round(deactivation, 2)}-psa{psa_press}.csv')

    # Save the results to a CSV file
    df = pd.DataFrame({key: arrays[key] for key in result_keys})

    output_dir = Path(__file__).resolve().parent.parent / "Results_Simulations"
    output_dir.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_dir / file_name, mode='w', index=False)

    return deactivation, cat_rate, tres

