###########################################################
# 2D MONTE CARLO NEUTRON TRANSPORT CODE #
###########################################################

import sys

sys.path.append("src")
from config.config import N_GENERATIONS
from src.reactor import Reactor

if __name__ == "__main__":

    reaction = Reactor(N_GENERATIONS)
    reaction.main_neutron_loop()
