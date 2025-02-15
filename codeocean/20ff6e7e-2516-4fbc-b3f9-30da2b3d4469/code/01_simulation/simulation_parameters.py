# Young moduli
E_hydrogel=2.1e4 # dECM fibirn hydrogel
E_PDMS=0.98e6 # From COMSOL library
E_cells=3e4

# Poisson coefficients
nu_cells=0.3 # Implies some water exchange
nu_hydrogel_self_condensation=0.2 # dECM fibrin hydrogel; under slow compression, water displacement is largely possible and the Poisson coefficient is expected to be low
nu_hydrogel_chip_compression=0.49 # These are rapid movements, so that the hydrogel is essentially incompressible. The 0.49 is for computational purposes
nu_PDMS=0.49 # From COMSOL library

# Gravity-related
rho=1000 # General density used. As we de not evaluate gravitational effects, this is without importance for now
g_acc=0 # No gravity effects

# Empirical parameters
isotropic_cell_contraction=0.3 # Where applicable, substantial cell self-contractility to enable for large potential self-condenstation effects
chip_compression=0.15 # Chip compression used for the experimental data acquisition


