Roberto S. Pinna 2024

# MC_2D_reactor
The Reactor class simulates the behavior of a nuclear reactor,
specifically the transport of neutrons within the reactor.
It uses a Monte Carlo method to track the interactions of neutrons
with the reactor's materials.

## Configuration file
User settings are in the config/config.py file, the configuration file for the nuclear reactor simulation codebase. It contains a set of variables that define the simulation parameters, which can be adjusted to change the behavior of the simulation.

Here's a breakdown of the variables defined in config.py:

N_GENERATIONS: The number of neutron generations to simulate.

ABS_COEFF: The boundary absorption coefficient, which represents the probability of a neutron being absorbed at the boundary of the reactor.

CHAIN: A boolean flag that determines whether the chain reaction is activated or not. If True, the simulation will continue to generate new neutrons until the specified number of generations is reached.

VERBOSE: A boolean flag that controls the verbosity of the simulation output. If True, the simulation will print out detailed information about the neutron interactions and reactor state.

LOG_DATA: A boolean flag that determines whether to log simulation data or not. If True, the simulation will write data to a log file.

TARGET_PWR: The target power flux for the PID (Proportional-Integral-Derivative) controller, which is used to regulate the reactor's power output.

PROP_K and INT_K: The proportional and integral gains for the PID controller, respectively.

These variables are used throughout the simulation code to control the behavior of the reactor and the simulation. By adjusting these values, users can change the simulation parameters to explore different scenarios or optimize the reactor's performance.

For example, increasing the value of N_GENERATIONS would simulate a longer period of time, while decreasing the value of ABS_COEFF would reduce the probability of neutron absorption at the boundary. Similarly, adjusting the TARGET_PWR value would change the desired power output of the reactor, which would affect the behavior of the PID controller.

## Running the Codebase

To run the codebase, follow these steps:

Clone the repository: Clone the codebase repository from GitHub or another version control system.
Install dependencies: Install the required dependencies listed in the requirements.txt file using pip:

pip install -r requirements.txt

Set up the configuration: 
Create a config.py file in the config/ directory of the codebase, or modify the existing one to set the desired simulation parameters.

Run the simulation: 
Run the simulation using the python command, specifying the main.py file as the entry point:

python reactor.py

This will start the simulation, and the code will begin to execute the neutron transport algorithm.
