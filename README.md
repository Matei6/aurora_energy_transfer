# Auroral Electron Energy Transfer Model

A student TIPE project that simulates electron collisions, excitation, and ionization in a simplified atmospheric model.

## Context

This project was completed in 2025 by Matei CRISTINAR as an individual project. The atmospheric composition source identified for the project is the MSIS-E90 model. The scripts use hard-coded atmospheric profiles and do not load MSIS-E90 data.

## Physical model

Altitude is expressed in kilometres and energy thresholds in electronvolts. Below 100 km, the atmospheric composition is set to 78% N2 and 22% O2. Above 100 km, the code uses these linear profiles, with proportions expressed in percent:

$$
x_{N_2}(h)=-0.38h+116,\qquad
x_{O_2}(h)=-0.14h+36,\qquad
x_O(h)=0.57h-57.
$$

The temperature profile is defined according to the US Standart atmosphere level

$$
T(h)=\begin{cases}
270.7-2.8(h-51), & 51\le h<71,\\
214.7-2.0(h-71), & 71\le h<85,\\
186.9, & 85\le h<91,\\
186.9+3.67(h-91), & 91\le h\le250.
\end{cases}
$$

Pressure is calculated by numerically integrating the coded hydrostatic integrand from the altitude to a layer reference altitude, then multiplying a hard-coded reference pressure by the exponential of the integral scaled by $10^3$. The integrand and gravity correction are:

$$
f(z)=\frac{l(z)gM(z)}{RT(z)},\qquad
l(h)=\left(\frac{6378}{6378+h}\right)^2.
$$

The mean molar mass $M(h)$ is calculated from the coded species proportions and molar-mass factors. The mean free path uses the temperature, pressure, and a composition-weighted radius:

$$
d(h)=10^{-2}\left(66\times10^{-12}x_O+146\times10^{-12}x_{O_2}+158\times10^{-12}x_{N_2}\right),\qquad
\lambda(h)=\frac{kT(h)10^{-3}}{\pi P(h)d(h)^2}.
$$

The coded excitation and ionization thresholds for O, O2, and N2 determine the energy subtracted from the electron at modeled collisions. The simulation does not load tabulated collision cross sections.

## Method

`src/aurora_tipe/averaged_ionization_simulation.py` advances an electron through altitude segments, calculates a Poisson mean $\mu=\Delta h/\lambda$, samples zero to four collisions per segment, places sampled collisions uniformly within the segment, selects an atmospheric species, and updates the electron energy at each collision. It computes an average ionization profile over 1,000 trajectories. The altitude and initial energy arguments are passed through to the trajectory simulations. A trajectory with no ionizations produces a zero-valued ionization profile.

## Repository structure

- `src/aurora_tipe/single_electron_simulation.py`: single-electron energy and altitude simulation.
- `src/aurora_tipe/excitation_ionization.py`: excitation and ionization helper functions.
- `src/aurora_tipe/electron_path_illustration.py`: randomized `turtle` path illustration.
- `src/aurora_tipe/atmospheric_model.py`: atmospheric profiles and mean-free-path plot.
- `src/aurora_tipe/averaged_ionization_simulation.py`: main averaged ionization simulation.
- `src/aurora_tipe/atmospheric_density_prototype.py`: atmospheric density and occupied-area experiment.
- `archive/collision_prototype.py`: archived collision prototype, preserved without modification.
- `figures/`: output directory for generated figures.
- `docs/35466_Cristinar_Matei_TIPE.pdf`: project PDF.
- `requirements.txt`: external Python dependencies.

## Installation and usage

From the `final_project` directory in Windows PowerShell:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python .\src\aurora_tipe\averaged_ionization_simulation.py --altitude 250 --energy 5000
```

`--altitude` is in kilometres and `--energy` is in electronvolts. Both options default to 250 km and 5,000 eV. The main script displays the ionization profile and saves it to `figures/ionization_profile.png`.

## Results

The main script generates an ionization-rate profile against altitude and saves the figure to `figures/ionization_profile.png`. It uses a logarithmic x-axis when the profile contains positive values and a linear x-axis for an all-zero profile. It does not save numerical data files.

## Limitations

- Atmospheric profiles and pressure reference values are hard-coded; the scripts do not load MSIS-E90 data.
- The collision model uses mean-free-path-based sampling and does not include tabulated species-dependent cross sections.
- The scripts contain duplicated model code and use different scaling factors in their mean-free-path expressions.
- The archived collision prototype is not part of the runnable simulation.
- Model values and approximations require scientific validation against the chosen data sources.

## References

MSIS-E90 atmospheric-composition model. US Standart atmosphere level.

## Author and license

Author: Matei CRISTINAR (TIPE, 2025).

License: MIT. See `LICENSE`.
