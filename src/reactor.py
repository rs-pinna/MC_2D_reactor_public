import sys

sys.path.append("src")

import time
import numpy as np
from scipy.stats import maxwell
from src.calculation import *
from src.utils import *
from src.constants import *
from config.config import *

np.random.seed()


class Reactor:
    """
    Class Overview
    The Reactor class simulates the behavior of a nuclear reactor,
    specifically the transport of neutrons within the reactor.
    It uses a Monte Carlo method to track the interactions of neutrons
    with the reactor's materials.
    - init: Initializes the reactor with a specified number of neutrons
    and sets various parameters, such as the chain reaction flag, absorption
    coefficient, and verbosity level.
    - main_neutron_loop: Simulates the neutron transport process, iterating over
    the specified number of neutrons and updating the reactor's state at each step.
    - neutron_life: Simulates the life of a single neutron, tracking its interactions
    with the reactor's materials and updating its position and energy.
    fuel, cladding, moderator: These methods simulate the behavior of the neutron
    within each of the reactor's regions (fuel, cladding, and moderator, respectively).
    - interact: Simulates the interaction of the neutron with the reactor's materials,
    including fission, capture, and inelastic scattering.
    """

    def __init__(self, Neutrons_Number: int) -> None:
        self.Neutrons_Number = Neutrons_Number
        self.chain = CHAIN
        self.abs_coeff = ABS_COEFF
        self.verbose = VERBOSE
        self.log_data = LOG_DATA
        self.alive = 1
        self.Leakage = 0
        self.max_d = 0
        self.real_d = 0
        self.Neutrons_Produced = 1
        self.Fission = 0
        self.Capture = 0
        self.Scattering = 0
        self.region = 0
        self.output_log = []
        self.decimals = int(-np.log10(X_STEP))
        self.active_pop = 100 / PITCH**2  # base neutron flux

    def main_neutron_loop(self):
        """
        Simulates the neutron transport loop in a nuclear reactor.

        This method performs the main simulation loop for neutron transport
        within the reactor. It iterates over a specified number of neutrons,
        updating their parameters and interactions within the reactor materials.

        The loop tracks the reactor's active neutron population and adjusts
        based on the effective multiplication factor (keff). It also calculates
        the probability of fission and power output, using a PID controller for
        power regulation.

        Neutron parameters such as energy, speed, position, and interaction
        type are initialized for each neutron, and the neutron life cycle is
        simulated.

        Attributes updated include:
        - start and end time of the simulation
        - reactor time
        - average delta time between generations
        - power output and energy output calculations

        Parameters
        ----------
        self : Reactor
            The instance of the Reactor class containing attributes related to
            the neutron transport simulation.

        Returns
        -------
        None
        """

        self.start_time = time.time()
        self.i = -1
        self.reactor_time = 0  # seconds
        iter_limit = 2 * self.Neutrons_Number + 30000
        while self.i < self.Neutrons_Number - 1:
            self.i += 1
            if self.i % 100 == 0:
                print(
                    f"Neutron {self.i} | {self.i/(time.time() - self.start_time + O_SMALL):.2f} n/s"
                )
            if self.i > iter_limit:
                # to avoid infinite loops
                print("Iter limit reached ", self.i)
                print("Neutrons projected ", self.Neutrons_Number)
                self.Neutrons_Number = self.i
                break

            # update of the active population by multiplication factor
            self.keff = (
                self.Neutrons_Produced + (1 - self.abs_coeff) * self.Leakage
            ) / (self.Capture + self.Fission + self.Leakage + O_SMALL)
            self.active_pop *= self.keff
            if self.active_pop <= 1 and CHAIN:
                print("Active population null")
                break

            # probability of fission and power output
            self.fission_prob = self.Fission / (
                self.Capture + self.Fission + self.Leakage * self.abs_coeff + O_SMALL
            )
            # power output per unit area
            self.energy_output = self.active_pop * self.fission_prob * URAN_FISSION_E
            self.delta_t = self.reactor_time / (
                self.i + 1
            )  # avg delta t between generations
            self.power_output = self.energy_output / (
                self.delta_t + O_SMALL
            )  # power flux of the current generation

            # power control
            pid_control(self)

            #### Initialize neutron parameters
            self.alive = 1  # Neutron life parameter
            self.interaction = 0  # Interaction control parameter
            self.intertype = 0  # Interaction type
            self.region = 0  # Region parameter [fuel,cladding,moderator]=[0,1,2]
            self.real_d = 0  # Travel distance parameter

            #### Initialize neutron location and energy
            # Maxwellian Neutron Energy Distribution [MeV]
            self.energy_i = maxwell.rvs()
            self.speed_i = 1.383e7 * np.sqrt(self.energy_i)  # neutron speed [m/s]
            # Random number for radial locations(powerlaw distribution)
            rnd = np.random.power(2)
            # Direction of neutron over 2pi
            self.r_i = rnd * R_FUEL
            self.theta_i = 2 * np.pi * np.random.random()
            # Initial coordinates of neutron location
            self.x_i = self.r_i * np.cos(self.theta_i)
            self.y_i = self.r_i * np.sin(self.theta_i)

            ###########
            # launch neutron life simulator
            self.neutron_life()
            ###########

        self.end_time = time.time()
        output_data(self)

    def neutron_life(self):
        """
        Simulates the neutron life cycle within the reactor.

        This method performs the main neutron life cycle loop. It iterates
        over a specified number of steps, updating the neutron's parameters
        and interactions within the reactor materials.

        The loop evaluates macroscopic cross sections, simulates the
        neutron's movement and interaction, and updates the neutron's
        parameters accordingly.

        Parameters
        ----------
        self : Reactor
            The instance of the Reactor class containing attributes related to
            the neutron transport simulation.

        Returns
        -------
        None
        """
        self.n_step = -1
        neutron_log(self)
        while self.alive == 1:
            # print(self.n_step)
            self.n_step += 1
            self.intertype = 0
            # evaluate macroscopic cross sections
            self.sig = CrossSections(self.energy_i, self.region)
            # maximum distance neutron goes before interaction
            self.max_d = np.round(
                np.abs((1 / self.sig[3]) * np.log(np.random.random())),
                decimals=self.decimals,
            )
            #### simulate regions
            if self.interaction == 1:
                self.interact()
                continue
            if self.region == 0:
                self.fuel()
                continue
            if self.region == 1:
                self.cladding()
                continue
            if self.region == 2:
                self.moderator()
                continue

    def fuel(self):
        """
        Simulates the neutron's interaction in the fuel region.

        This method calculates the macroscopic cross sections in the fuel
        region and simulates the neutron's movement and interaction. If the
        neutron is leaving the fuel region, it updates the neuron's parameters
        accordingly.

        Parameters
        ----------
        self : Reactor
            The instance of the Reactor class containing attributes related to
            the neutron transport simulation.

        Returns
        -------
        None
        """
        self.A = URAN_A  # Mass Number(A) of Uranium
        # self.density = URAN_RHO  # g / cc

        # Check distance for nearest surface
        # ------------------------------------------
        df, dmin = d_fuel_surfaces(self)

        if self.max_d >= dmin:
            # Neutron is leaving the fuel region(0)
            # and entering the cladding region(1)
            if self.verbose:
                print("fuel -> cladding")
            self.real_d = dmin + 4 * X_STEP
            self.intertype = 4
            self.region = 1
        else:
            # Neutron is interacting at fuel region
            self.real_d = self.max_d
            self.interaction = 1

        update_position(self, self.x_i, self.y_i, self.real_d, self.theta_i)

        return None

    def cladding(self):
        """
        Simulates the neutron's interaction in the cladding region.

        This method calculates the macroscopic cross sections in the cladding
        region and simulates the neutron's movement and interaction. If the
        neutron is leaving the cladding region, it updates the neuron's
        parameters accordingly.

        Parameters
        ----------
        self : Reactor
            The instance of the Reactor class containing attributes related to
            the neutron transport simulation.

        Returns
        -------
        None
        """
        # =========================================================#
        # Cladding Region #
        # =========================================================#
        self.A = ZIRC_A  # Mass Number(A) of Zirconium
        # self.density = ZIRC_RHO  # g / cc

        # Check distance for nearest surface

        dc, dmin = d_cladding_surfaces(self)

        if self.max_d >= dmin:
            # Neutron is leaving the cladding region(1)
            self.intertype = 4
            self.real_d = dmin + 4 * X_STEP
            if dmin == dc[0] or dmin == dc[1]:
                if self.verbose:
                    print("cladding -> fuel")
                self.region = 0
            elif dmin == dc[2] or dmin == dc[3]:
                if self.verbose:
                    print("cladding -> moderator")
                self.region = 2
        else:
            self.real_d = self.max_d
            # Neutron is interacting at cladding region
            if self.verbose:
                print("interaction in cladding")
            self.interaction = 1

        update_position(self, self.x_i, self.y_i, self.real_d, self.theta_i)

        return None

    def moderator(self):
        """
        Simulates the neutron's interaction in the moderator region.

        This method calculates the macroscopic cross sections in the moderator
        region and simulates the neutron's movement and interaction. If the
        neutron is leaving the moderator region, it updates the neutron's
        parameters accordingly.

        Parameters
        ----------
        self : Reactor
            The instance of the Reactor class containing attributes related to
            the neutron transport simulation.

        Returns
        -------
        None
        """
        self.A = WATER_A  # Mass Number(A) of H
        # self.density = WATER_RHO  # g/cc

        dm, dmin = d_moderator_surfaces(self)

        if self.max_d >= dmin:
            # Neutron is leaving the moderator region(2)
            self.real_d = dmin
            self.intertype = 4
            update_position(self, self.x_i, self.y_i, self.real_d, self.theta_i)

            if dmin == dm[0] or dmin == dm[1]:
                if self.verbose:
                    print("moderator -> cladding")
                self.region = 1
            elif dmin == dm[2] or dmin == dm[4]:
                if self.verbose:
                    print("leaving from left/right boundary")
                periodic_boundary(self, on="x")

            elif dmin == dm[3] or dmin == dm[5]:
                if self.verbose:
                    print("leaving from top/bottom boundary")
                periodic_boundary(self, on="y")
        else:
            self.real_d = self.max_d
            if self.verbose:
                print("interaction in moderator")
            self.interaction = 1
            update_position(self, self.x_i, self.y_i, self.real_d, self.theta_i)

        return None

    def interact(self):
        """
        Handles the interaction of a neutron.

        This method randomly determines the interaction type based on the
        macroscopic cross sections of the material. The interaction type can
        be fission, absorption, or inelastic scattering.

        If the interaction type is fission, the function calls `fission`.
        If the interaction type is absorption, the function calls `capture`.
        If the interaction type is inelastic scattering, the function calls
        `inelastic_scattering`.

        The function then logs the neutron's data if logging is enabled.

        Parameters
        ----------
        self : Reactor
            The instance of the Reactor class containing attributes related to
            the neutron transport simulation.

        Returns
        -------
        None
        """

        self.real_d = 0  # neutron interact on the same spot
        self.interaction = 0
        rnd = np.random.random()
        self.reactor_time += O_SMALL

        if rnd <= (self.sig[0] / self.sig[3]):
            fission(self)
        elif (self.sig[0] / self.sig[3]) < rnd and (
            rnd <= (self.sig[0] + self.sig[1]) / self.sig[3]
        ):
            capture(self)
        else:
            inelastic_scattering(self)
        #
        if self.log_data:
            neutron_log(self)

        return None
