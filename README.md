# Catalytic Methane Pyrolysis: Process Simulation and Techno-Economic Assessment

This repository contains the code, process simulations and results for the research study:

> **Techno-Economic Comparison of Ni-, Fe-, and Activated-Carbon-Based Catalysts for Methane Pyrolysis in a Fluidized-Bed Reactor**
>
> Gabriele Patruno\*, Panji B. Tamarona\*, Thijs J. H. Vlugt, and Mahinder Ramdin
>
> \*These authors contributed equally to this work.

It compares hydrogen production by catalytic methane pyrolysis (CH₄ → C + 2 H₂) on three catalysts:

| Folder | Catalyst | Base-case temperature | Base-case catalyst inventory |
|--------|----------|-----------------------|------------------------------|
| [`catalysts/ni/`](catalysts/ni) | Nickel (Ni) | 600 °C | 950 kg |
| [`catalysts/fe/`](catalysts/fe) | Iron (Fe) | 725 °C | 10 000 kg |
| [`catalysts/ac/`](catalysts/ac) | Activated carbon (AC) | 825 °C | 210 000 kg |

The pyrolysis plant produces 100 t/d of H₂. It has a fluidised-bed reactor (modelled as a CSTR or a PFR with catalyst circulation and deactivation), cyclone carbon separation, compression, and PSA purification. Each catalyst is evaluated with four reactor heating strategies:

- methane combustion (`ch4`)
- methane combustion with CO₂ capture and storage (CCS)
- hydrogen combustion (`h2`)
- electric heating (`electric`)

Steam methane reforming (SMR) and water electrolysis are compared using literature values, so they have no models in this repository.

## Repository structure

Every catalyst folder has the same layout:

```
catalysts/
└── ni/ | fe/ | ac/
    ├── aspen/                    # Aspen Plus flowsheets: {cstr,pfr}-{ch4,h2,electric}-<Cat>.apw
    ├── simulation/
    │   ├── kinetic_models.py     # Isothermal CSTR and PFR reactor models (kinetics + deactivation)
    │   ├── aspen_simulations.py  # Aspen Plus COM interface: couples the reactor model to the flowsheet
    │   ├── run_simulations.ipynb # Simulation configuration and entry point
    │   └── keq_aspen.py          # (Fe only) fit of the equilibrium constant to the Aspen Kp form
    ├── tea/                      # Techno-economic assessment notebooks (OpenPyTEA)
    │   ├── CH4_CCS.ipynb         # CH4 heating with CCS: base case, sensitivities, Monte Carlo
    │   ├── CH4.ipynb             # CH4 heating without CCS
    │   ├── H2.ipynb              # H2 heating
    │   └── Electric.ipynb        # Electric heating
    └── results/
        ├── simulations/          # Process simulation results, one CSV per case
        └── tea/                  # TEA results (JSON)
```

### Reactor models

The `kinetic_models.py` files implement isothermal **CSTR** and **PFR** models of the pyrolysis reactor. Each model includes:

- the catalyst-specific rate law (Arrhenius / Langmuir-Hinshelwood type, depending on the catalyst)
- time-dependent catalyst deactivation, integrated with RK4 over the catalyst residence time
- an equilibrium limit on conversion

Each simulation solves for the catalyst circulation rate that gives the target average catalyst deactivation (`mode: constant_da`). The resulting reactor conversion is written into the Aspen Plus flowsheet. The flowsheet is then rerun until the recycle loop converges.

The kinetic and deactivation parameters are set in `sim_config` in each `run_simulations.ipynb`.

### Simulation results

Each file in `results/simulations/` is named

```
sim-<CSTR|PFR>-<capacity>TPD-<heating>-t<T °C>-p<P bar>-c<catalyst kg>-da<deactivation>-psa<PSA bar>.csv
```

and contains the stream flows, catalyst rate and residence time, reactor size and count, compressor/turbine powers, heat-exchanger and cooler areas, utility and fuel use, and CO₂ emissions. The TEA notebooks use these values to size and cost the equipment.

### TEA results

| File in `results/tea/` | Content |
|------------------------|---------|
| `LCOH_heatings.json` | Levelised cost of hydrogen (LCOH) per heating strategy |
| `CO2_results.json` | Specific CO₂ emissions per heating strategy |
| `LCOH_temperature.json` | LCOH vs. reactor temperature (CH₄ + CCS) |
| `LCOH_avg.json` | Levelised cost breakdown |
| `Direct_costs_avg.json`, `Fixed_capital_avg.json` | Equipment and fixed-capital cost breakdown |
| `Variable_opex_avg.json`, `Fixed_opex_avg.json` | Operating-cost breakdown |

The Monte Carlo results (`LCOH_MonteCarlo.json`) are not included. Each file is 130–200 MB, and GitHub rejects files over 100 MB. To regenerate them, uncomment the Monte Carlo cells at the end of `tea/CH4_CCS.ipynb`. These runs are slow.

## Where each figure's data comes from

The repository has no plotting scripts. This table shows which results each figure in the study is based on.

| Figure | Content | Data |
|--------|---------|------|
| 3 | Baseline comparison of the three catalysts | `catalysts/*/results/tea/*_avg.json` |
| 4 | Catalyst circulation rate | `catalysts/*/results/simulations/` (base cases, `Cat Rate [kg/h]`) |
| 5 | LCOH vs. reactor temperature | `catalysts/*/results/tea/LCOH_temperature.json` |
| 6 | Monte Carlo LCOH distributions | `LCOH_MonteCarlo.json` (regenerate, see above) |
| 7 | [TODO] | [TODO] |
| 8 | LCOH and CO₂ emissions per heating strategy | `catalysts/*/results/tea/LCOH_heatings.json`, `CO2_results.json` |
| S1 | AC kinetic parameter fit | parameters in `catalysts/ac/simulation/run_simulations.ipynb` |
| S2 | [TODO] | [TODO] |
| S3–S5 | LCOH, CSTR vs. PFR (Ni, Fe, AC) | `catalysts/{ni,fe,ac}/tea/CH4_CCS.ipynb` |
| S6–S8 | Cost breakdown (Ni, Fe, AC) | `catalysts/{ni,fe,ac}/results/tea/*_avg.json` |
| S9, S10 | Fixed capital investment, fixed OPEX | `catalysts/*/results/tea/Fixed_capital_avg.json`, `Fixed_opex_avg.json` |
| S11–S13 | Temperature sensitivity (Ni, Fe, AC) | `catalysts/{ni,fe,ac}/results/tea/LCOH_temperature.json` and the temperature-sweep CSVs |
| S14–S16 | Tornado sensitivity (Ni, Fe, AC) | computed in `catalysts/{ni,fe,ac}/tea/CH4_CCS.ipynb` |

## Installation

Python ≥ 3.12 is required (developed with Python 3.14).

```bash
git clone https://github.com/pbtamarona/CatalyticMethanePyrolysis.git
cd CatalyticMethanePyrolysis
python -m venv .venv
.venv\Scripts\activate          # Windows  (source .venv/bin/activate on Linux/macOS)
pip install .
```

This installs the dependencies listed in `pyproject.toml`, including [OpenPyTEA](https://github.com/pbtamarona/OpenPyTEA) at the commit used for the study. To reproduce the exact environment, run `pip install -r requirements.txt` instead.

**Process simulations** need Windows and a licensed copy of **Aspen Plus**. The code controls Aspen Plus through COM, using `pywin32`. The **TEA notebooks** run on any platform without Aspen Plus, because they read the CSVs that are already in `results/simulations/`.

## Usage

### 1. Run a process simulation (optional; requires Aspen Plus)

Open `catalysts/<cat>/simulation/run_simulations.ipynb` and edit `sim_config` for the case you want:

```python
sim_config = {
    'plant_capacity': 100,                       # t/d H2
    'reactor': {'type': 'CSTR',                  # 'CSTR' or 'PFR'
                'temperature': 600,              # °C
                'pressure': 1.07,                # bar
                'heating': 'ch4'},               # 'ch4', 'h2' or 'electric'
    'catalyst': {'weight': 950, 'density': 850}, # kg, kg/m3
    'mode': {'type': 'constant_da', 'value': 0.2, ...},
    ...
}
run_scenario(sim_config)
```

The code finds the flowsheet in `../aspen/` and writes the result CSV to `../results/simulations/`. Both paths are relative to `aspen_simulations.py`, so the simulation works no matter which folder you start it from.

### 2. Run the techno-economic assessment

For each catalyst, run the notebooks in `catalysts/<cat>/tea/` in this order:

1. `CH4_CCS.ipynb` writes the base case, cost breakdowns, temperature and sensitivity results, and creates `LCOH_heatings.json` and `CO2_results.json`.
2. `CH4.ipynb`, `H2.ipynb` and `Electric.ipynb` add their heating strategy to `LCOH_heatings.json` and `CO2_results.json`.

## License

- **Code** (`*.py`, `*.ipynb`): [MIT](LICENSE)
- **Data** (`results/`, Aspen Plus files in `aspen/`): [CC-BY-4.0](LICENSE-DATA)
