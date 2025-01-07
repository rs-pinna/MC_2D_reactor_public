from src.constants import *
import pandas as pd
import matplotlib.pyplot as plt


def output_data(self: object) -> None:
    """
    Processes and logs neutron transport data.

    This function prints output data and, if logging is enabled, logs the
    current state of the neutron transport simulation. It stores the data
    in a DataFrame, removes duplicates, and exports the log to a CSV file.

    Parameters
    ----------
    self : object
        The instance of the class containing attributes related to the
        neutron transport simulation.

    Returns
    -------
    None
    """

    output_print(self)
    if self.log_data:
        neutron_log(self)
        self.output_log = pd.DataFrame(
            self.output_log,
            columns=[
                "generation",
                "step",
                "time",
                "region",
                "x",
                "y",
                "theta",
                "real_d",
                "speed",
                "energy",
                "intertype",
                "active_pop",
                "power_output",
                "abs_coeff",
                "keff",
                "fission_prob",
            ],
        )
        self.output_log.drop_duplicates(inplace=True)
        self.output_log.to_csv("data/output_log.csv", float_format="%.6e", index=False)
        plotting(self)


def neutron_log(self: object) -> None:

    # append data to list output_log
    """
    Logs the current state of the neutron transport simulation.

    This method records the current attributes of the neutron transport
    simulation into the `output_log` list. These attributes include the
    generation index, step count, reactor time, region, neutron position
    (x, y), direction angle, real distance traveled, speed, energy,
    interaction type, active population, power output, and absorption
    coefficient.

    Parameters
    ----------
    self : object
        The instance of the class containing attributes related to the
        neutron transport simulation.

    Returns
    -------
    None
    """

    values = [
        self.i,
        self.n_step,
        self.reactor_time,
        self.region,
        self.x_i,
        self.y_i,
        self.theta_i,
        self.real_d,
        self.speed_i,
        self.energy_i,
        self.intertype,
        self.active_pop,
        self.power_output,
        self.abs_coeff,
        self.keff,
        self.fission_prob,
    ]

    self.output_log.append(values)

    return None


def output_print(self: object) -> None:
    """
    Prints the results of the neutron transport simulation.

    This function prints the results of the neutron transport simulation to
    the console. It calculates and displays the number of simulated
    generations, calculation speed, reactor time, average neutron duration,
    number of interactions, scattering events, absorption events, capture
    events, fission events, average nu, number of neutrons produced by fission,
    number of neutrons leaked from the system, number of neutrons leaked into
    the system, effective multiplication factor, and average fission
    probability.

    Parameters
    ----------
    self : object
        The instance of the class containing attributes related to the
        neutron transport simulation.

    Returns
    -------
    None
    """
    Absorption = self.Fission + self.Capture
    Leakage_in = (
        1 - self.abs_coeff
    ) * self.Leakage  # Leakage  # neutrons leaked into the system
    nu = self.Neutrons_Produced / (self.Fission + O_SMALL)
    Interactions = self.Scattering + Absorption

    print(
        "Number of stories.......................= ", self.Neutrons_Number
    )  # Number of simulated generations
    print(
        f"Calculation speed.(neutrons/s)...........= {self.Neutrons_Number/(self.end_time - self.start_time + O_SMALL):.1f}"
    )
    print(f"Reactor time (seconds)...................= {self.reactor_time:.2f}")
    print(
        f"Avg neutron duration (seconds)...........= {self.reactor_time/self.Neutrons_Number:.3f}"
    )
    print(f"Number of Interactions...................= {Interactions}")
    print(f"Number of Scattering Events..............= {self.Scattering}")
    print(f"Number of Absorption Events..............= {Absorption}")
    print(f"  Number of Capture Events...............= {self.Capture}")
    print(f"  Number of Fission Events...............= {self.Fission}")
    print(
        f"Average nu...............................= {nu:.3f}"
    )  # Average number of neutrons produced per fission
    print(f"Number of Neutrons Produced by Fission...= {self.Neutrons_Produced}")
    print(f"Number of Neutrons Leaked from System....= {self.Leakage}")
    print(f"Number of Neutrons Leaked into System....= {Leakage_in:.0f}")
    print(f"Effective Multiplication Factor(keff)....= {self.keff:.6f}")
    print(f"Avg fission probability..................= {self.fission_prob:.3f}")


def plotting(self: object) -> None:
    """
    Plot the neutron trajectories and the power output of the reactor.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    The function first removes the duplicate points from the log of the simulation
    and then plots the neutron trajectories using a line plot. The colors of the
    lines are determined by the energy of the neutron at each point. The function
    also plots the power output of the reactor over time in a second figure.
    """

    plot_log = self.output_log.copy()
    plot_log.drop_duplicates(subset=["x", "y"], inplace=True)

    plot_log = plot_log.iloc[0:500]  # limit the number of lines in the plot
    plt.figure(figsize=(10, 10))
    plt.gcf().gca().add_artist(plt.Circle((0, 0), R_FUEL, fill=False))
    plt.gcf().gca().add_artist(plt.Circle((0, 0), R_CLAD_IN, fill=False))
    plt.gcf().gca().add_artist(plt.Circle((0, 0), R_CLAD_OUT, fill=False))
    plt.gcf().gca().add_artist(
        plt.Rectangle((-PITCH / 2, -PITCH / 2), width=PITCH, height=PITCH, fill=False)
    )

    color_map = plt.get_cmap("jet")
    for i in range(plot_log.shape[0] - 1):
        if (plot_log.iloc[i].x != plot_log.iloc[i + 1].x) and (
            plot_log.iloc[i].y != plot_log.iloc[i + 1].y
        ):
            plt.plot(
                [plot_log.iloc[i].x, plot_log.iloc[i + 1].x],
                [plot_log.iloc[i].y, plot_log.iloc[i + 1].y],
                c=color_map(plot_log.iloc[i].energy),
                alpha=0.5,
            )

    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.grid()
    plt.show()

    ##############
    plt.figure()
    plt.plot(self.output_log["time"], self.output_log["power_output"])
    plt.yscale("log")
    plt.xlabel("Time [s]")
    plt.ylabel("Power flux [W/m^2]")
    plt.grid()
    plt.show()
    return None
