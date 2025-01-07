import numpy as np

# =========================================================#
# Limits #
# =========================================================#

X_STEP = 1e-6
O_BIG = 1e10
O_SMALL = 1e-10
THERMAL_E = 2.5e-8  # energy of thermalization [MeV]

# =========================================================#
# Materials #
# =========================================================#
URAN_A = 238.02891  # Mass Number(A) of Uranium
URAN_RHO = 10.97  # g / cc
URAN_FISSION_E = 3.2e-11  # avg energy per fission [J]

ZIRC_A = 26.981539  # Mass Number(A) of Zirconium
ZIRC_RHO = 2.70  # g / cc

WATER_A = 1.00794  # Mass Number(A) of H
WATER_RHO = 1  # g/cc

# =========================================================#
# Geometry #
# =========================================================#
R_FUEL = 0.1  # Fuel radius [m]
R_CLAD_IN = 0.1  # Cladding inner radius [m]
R_CLAD_OUT = 0.15  # Cladding outer radius [m]
T_CLAD = R_CLAD_OUT - R_CLAD_IN  # Cladding thickness [m]
PITCH = 2  # Cell side length [m]

# =========================================================#
# Neutron energy groups #
# =========================================================#
# Unit is MeV
ENERGY_GROUPS = np.array(
    [
        3.0e-08,
        3.0e-07,
        3.0e-06,
        3.0e-05,
        3.0e-04,
        3.0e-03,
        3.0e-02,
        3.0e-01,
        3.0e00,
        3.0e01,
    ]
)

# =========================================================#
# Energy Dependent Macroscopic Cross Section Data #
# =========================================================#
# Unit is 1/cm
# sigma_x[EnergyGroup][Region]
# coluns are regions: fuel, cladding, moderator

# fission cross-section
SIGMA_F = np.array(
    [
        [25.8, 0.0, 0.0],
        [9.3, 0.0, 0.0],
        [1.28, 0.0, 0.0],
        [0.107, 0.0, 0.0],
        [0.25, 0.0, 0.0],
        [0.246, 0.0, 0.0],
        [0.106, 0.0, 0.0],
        [0.0602, 0.0, 0.0],
        [0.0596, 0.0, 0.0],
        [0.105, 0.0, 0.0],
    ]
)

# capture cross-section
SIGMA_C = np.array(
    [
        [4.30e00, 3.95e-04, 5.76e-06],
        [2.12e00, 1.26e-04, 1.84e-06],
        [2.51e-01, 3.98e-05, 6.02e-07],
        [9.90e-02, 1.25e-05, 2.02e-07],
        [4.28e-02, 4.39e-06, 1.27e-07],
        [8.23e-02, 1.06e-05, 2.24e-07],
        [3.29e-02, 4.52e-06, 6.63e-07],
        [1.10e-02, 2.83e-04, 2.56e-07],
        [1.34e-03, 7.83e-03, 3.34e-06],
        [1.41e-06, 1.71e-02, 3.34e-06],
    ]
)

# scattering cross-section
SIGMA_S = np.array(
    [
        [4.36e01, 2.41e-01, 6.91e-01],
        [9.79e00, 2.40e-01, 6.86e-01],
        [2.45e00, 2.31e-01, 6.83e-01],
        [2.30e00, 2.23e-01, 6.82e-01],
        [1.52e00, 2.14e-01, 6.81e-01],
        [9.38e-01, 2.06e-01, 6.69e-01],
        [6.88e-01, 2.66e-01, 5.72e-01],
        [4.77e-01, 3.44e-01, 2.65e-01],
        [3.88e-01, 1.76e-01, 7.36e-02],
        [2.76e-01, 1.44e-01, 1.27e-02],
    ]
)
