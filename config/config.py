# =========================================================#
# User Config #
# =========================================================#
N_GENERATIONS = 1000  # number of neutron generations
ABS_COEFF = 0.5  # boundary absorption coefficient (alias control rod position)
CHAIN = False  # activate the chain reaction (the loop may never finish)
VERBOSE = False  # activate printout
LOG_DATA = True  # activate logging
TARGET_PWR = 1e6  # target power flux for the PID [W/m^-2]
PROP_K = 1e-2  # pid coefficient
INT_K = 1e-2  # pid coefficient
