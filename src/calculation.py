import sys

sys.path.append("src")
from src.constants import *
from config.config import *
from src.utils import neutron_log
import numpy as np


def CrossSections(Energy: float, region: int) -> list:
    """
    Calculates the macroscopic cross sections for a given energy and region.

    Parameters
    ----------
    Energy : float
        The energy of the neutron in MeV.
    region : int
        The region of the reactor.

    Returns
    -------
    list
        A list containing the macroscopic cross sections for fission, capture, scattering, and total.
    """
    sig_f = np.interp(Energy, ENERGY_GROUPS, SIGMA_F[:, region])
    sig_c = np.interp(Energy, ENERGY_GROUPS, SIGMA_C[:, region])
    sig_s = np.interp(Energy, ENERGY_GROUPS, SIGMA_S[:, region])
    sig_t = sig_f + sig_c + sig_s
    sig = [sig_f, sig_c, sig_s, sig_t]

    return sig


def d_fuel_surfaces(self: object) -> list:
    """
    Calculate the distances to the fuel surface from a point (x, y) with direction angle theta.

    Parameters
    ----------
    self : object
        Object containing the neutron's parameters

    Returns
    -------
    df : list
        Distances to the fuel surface
    dmin : float
        Minimum distance to the fuel surface

    Notes
    -----
    The grazing control is implemented to avoid singularities in the calculation of the distance
    to the fuel surface. If the neutron is too close to the fuel surface, it is moved slightly away
    from the surface to avoid the singularity. The distance to the fuel surface is then calculated
    taking into account this correction.

    The function returns a list of distances to the fuel surface and the minimum distance to the
    fuel surface.

    Raises
    ------
    ValueError
        If the neutron is outside the fuel region, or if the calculation of the distance to the
        fuel surface fails.
    """

    #############
    # Grazing control
    r_i = np.sqrt(self.x_i**2 + self.y_i**2)
    if np.abs(R_FUEL - r_i) < 8 * X_STEP:
        if self.verbose:
            print("Grazing impact on fuel")
        # Ensure that the neutron is inside the region
        self.x_i -= 8 * X_STEP * np.sign(self.x_i)
        self.y_i -= 8 * X_STEP * np.sign(self.y_i)
        self.x_i = np.round(self.x_i, decimals=self.decimals)
        self.y_i = np.round(self.y_i, decimals=self.decimals)
    #############

    a = 1
    b = 2 * (self.x_i * np.cos(self.theta_i) + self.y_i * np.sin(self.theta_i))
    c = self.x_i**2 + self.y_i**2 - R_FUEL**2
    delta = b**2 - 4 * a * c

    if delta < 0:
        raise ValueError("delta cannot be negative")

    sqrt_delta = delta**0.5
    df = [(-b - sqrt_delta) / (2 * a), (-b + sqrt_delta) / (2 * a)]

    df = np.round(df, decimals=self.decimals)
    d_pos = [i for i in (df) if i > O_SMALL]
    if not d_pos:
        raise ValueError("No valid positive distances")

    dmin = min(d_pos)

    if dmin > 2.1 * R_FUEL:
        raise ValueError("wrong distance")

    return df, dmin


def d_cladding_surfaces(self: object) -> list:
    """
    Calculates the distances to the inner and outer cladding surfaces
    from the current neutron position (self.x_i, self.y_i) and direction
    angle (self.theta_i).

    The function first ensures that the neutron is not grazing the cladding
    surfaces, adjusting its position if necessary. It then calculates the
    roots of the quadratic equations representing intersections with the
    inner and outer cladding surfaces. The distances are computed for both
    intersections and returned along with the minimum valid distance.

    Returns:
        dc (list): A list containing the computed distances to the inner
                   and outer cladding surfaces. Indices 0 and 1 correspond
                   to the inner cladding, and indices 2 and 3 correspond
                   to the outer cladding.
        dmin (float): The minimum valid positive distance to a cladding
                      surface.

    Raises:
        ValueError: If delta_out is negative or if the minimum distance
                    exceeds the maximum allowable distance.
    """

    #######
    # Grazing control
    r_i = np.sqrt(self.x_i**2 + self.y_i**2)
    if np.abs(r_i - R_CLAD_OUT) < 4 * X_STEP:
        if self.verbose:
            print("Grazing impact on outer cladding")
        # Ensure that the neutron is inside the region
        self.x_i -= 4 * X_STEP * np.sign(self.x_i)
        self.y_i -= 4 * X_STEP * np.sign(self.y_i)
        self.x_i = np.round(self.x_i, decimals=self.decimals)
        self.y_i = np.round(self.y_i, decimals=self.decimals)
    elif np.abs(r_i - R_CLAD_IN) < 4 * X_STEP:
        if self.verbose:
            print("Grazing impact on inner cladding")
        # Ensure that the neutron is inside the region
        self.x_i += 2 * X_STEP * np.sign(self.x_i)
        self.y_i += 2 * X_STEP * np.sign(self.y_i)
        self.x_i = np.round(self.x_i, decimals=self.decimals)
        self.y_i = np.round(self.y_i, decimals=self.decimals)
    ###########
    a = 1
    b = 2 * (self.x_i * np.cos(self.theta_i) + self.y_i * np.sin(self.theta_i))
    c_in = self.x_i**2 + self.y_i**2 - R_CLAD_IN**2
    c_out = self.x_i**2 + self.y_i**2 - R_CLAD_OUT**2

    delta_in = b**2 - 4 * a * c_in
    delta_out = b**2 - 4 * a * c_out
    dc = [0] * 4
    d_pos = []

    if delta_in >= 0:
        dc[0] = (-b - delta_in**0.5) / (2 * a)
        dc[1] = (-b + delta_in**0.5) / (2 * a)
    else:
        dc[0] = O_BIG
        dc[1] = O_BIG
    #
    if delta_out >= 0:
        dc[2] = (-b - delta_out**0.5) / (2 * a)
        dc[3] = (-b + delta_out**0.5) / (2 * a)
    else:
        raise ValueError("delta_out cannot be negative")

    dc = np.round(dc, decimals=self.decimals)
    d_pos = [d for d in (dc) if d > O_SMALL]
    dmin = min(d_pos)

    if dmin > 2.001 * R_CLAD_OUT:
        raise ValueError("wrong distance")

    return dc, dmin


def d_moderator_surfaces(self: object) -> list:
    """
    Calculates the distances to the boundaries of the moderator region
    from the current neutron position (self.x_i, self.y_i) and direction
    angle (self.theta_i).

    The function includes a grazing control mechanism to ensure the neutron
    is not too close to the outer cladding surface. It adjusts the neutron's
    position if necessary and calculates the distances to potential intersection
    points with the moderator boundaries and outer cladding surface.

    Returns:
        dm (list): A list containing the computed distances to various
                   boundaries. Indices 0 and 1 correspond to the outer cladding
                   surface, while indices 2 to 5 correspond to the right, top,
                   left, and bottom boundaries, respectively.
        dmin (float): The minimum valid positive distance to a boundary.

    Raises:
        ValueError: If the minimum distance exceeds the maximum allowable
                    distance.
    """
    #######
    # Grazing control
    r_i = np.sqrt(self.x_i**2 + self.y_i**2)
    if np.abs(r_i - R_CLAD_OUT) < 4 * X_STEP:
        if self.verbose:
            print("Grazing impact on moderator")
        self.x_i += 4 * X_STEP * np.sign(self.x_i)
        self.y_i += 4 * X_STEP * np.sign(self.y_i)
        self.x_i = np.round(self.x_i, decimals=self.decimals)
        self.y_i = np.round(self.y_i, decimals=self.decimals)
    #
    #########

    a = 1
    b = 2 * (self.x_i * np.cos(self.theta_i) + self.y_i * np.sin(self.theta_i))
    c_out = self.x_i**2 + self.y_i**2 - R_CLAD_OUT**2
    delta_out = b**2 - 4 * a * c_out
    dm = [0] * 6
    d_pos = []
    if delta_out >= 0:
        dm[0] = (-b - delta_out**0.5) / (2 * a)
        dm[1] = (-b + delta_out**0.5) / (2 * a)
    else:
        dm[0] = O_BIG
        dm[1] = O_BIG

    dm[2] = (PITCH / 2 - self.x_i) / np.cos(self.theta_i)  #
    # distance to right boundary
    dm[3] = (PITCH / 2 - self.y_i) / np.sin(self.theta_i)  #
    # distance to top boundary
    dm[4] = (-PITCH / 2 - self.x_i) / np.cos(self.theta_i)  #
    # distance to left boundary
    dm[5] = (-PITCH / 2 - self.y_i) / np.sin(self.theta_i)  #
    # distance to bottom boundary

    dm = np.round(dm, decimals=self.decimals)
    d_pos = [i for i in (dm) if i > O_SMALL]
    dmin = min(d_pos)

    if dmin > 1.001 * np.sqrt(2) * PITCH:
        raise ValueError("wrong distance")

    return dm, dmin


def update_position(self, x_old: float, y_old: float, d: float, theta_i: float) -> None:
    """
    Updates the neutron's position and time based on a given distance and direction angle.

    Parameters
    ----------
    x_old : float
        The neutron's old x-coordinate.
    y_old : float
        The neutron's old y-coordinate.
    d : float
        The distance to move the neutron.
    theta_i : float
        The direction angle of the neutron's movement.

    Notes
    -----
    The function includes a boundary check to ensure the neutron is not moved outside the
    reactor's boundaries. It also logs the neutron's movement if self.log_data is True.
    """
    x_new = np.round(x_old + d * np.cos(theta_i), decimals=self.decimals)
    y_new = np.round(y_old + d * np.sin(theta_i), decimals=self.decimals)

    ########
    # boundary check
    limit = PITCH / 2 - 4 * X_STEP
    x_new = np.clip(x_new, -limit, limit)
    y_new = np.clip(y_new, -limit, limit)
    ########

    ### update position and time
    self.reactor_time += d / self.speed_i
    self.r_i = np.sqrt(x_new**2 + y_new**2)
    self.x_i = np.round(x_new, decimals=self.decimals)
    self.y_i = np.round(y_new, decimals=self.decimals)
    if self.log_data:
        neutron_log(self)

    return None


def periodic_boundary(self: object, on: str) -> None:
    """
    Handles neutron boundary crossing in a periodic manner.

    Parameters
    ----------
    on : str
        The axis along which the neutron has crossed the boundary ('x' or 'y').

    Notes
    -----
    The function randomly determines if a neutron is absorbed based on the absorption
    coefficient. If not absorbed, the neutron's position is updated to reflect its
    reappearance on the opposite boundary. The function also updates the reactor time
    and logs the neutron's data if logging is enabled.

    Raises
    ------
    NotImplementedError
        If the specified boundary axis is not 'x' or 'y'.
    """

    self.alive = np.random.choice([0, 1], p=[self.abs_coeff, (1 - self.abs_coeff)])
    self.Leakage += 1

    ### update position and time
    self.reactor_time += 4 * X_STEP / self.speed_i
    if on == "x":
        self.x_i = np.round(-self.x_i * (1 - 4 * X_STEP), decimals=self.decimals)
    elif on == "y":
        self.y_i = np.round(-self.y_i * (1 - 4 * X_STEP), decimals=self.decimals)
    else:
        raise NotImplementedError("Dimension not implemented")

    if self.log_data:
        neutron_log(self)


def fission(self: object) -> None:
    """
    Handles the fission interaction of a neutron.

    Notes
    -----
    This function updates the interaction type to fission, sets the neutron as dead,
    increments the fission count, and determines the number of neutrons produced in
    the fission event. If chain reactions are enabled, it also adjusts the total
    number of neutrons accordingly.

    The function sets the color to red to represent fission events.
    """

    self.intertype = 1  # fission
    self._color = "r"
    self.Fission += 1
    self.alive = 0
    self.nf = np.random.choice([2, 3])
    self.Neutrons_Produced += self.nf
    ######### chain reaction  ########
    if self.chain:
        self.Neutrons_Number += self.nf
    ###################################


def capture(self: object) -> None:
    """
    Handles the capture interaction of a neutron.

    Notes
    -----
    This function updates the interaction type to capture, sets the neutron as dead,
    and increments the capture count.

    The function sets the color to blue to represent capture events.
    """
    self.intertype = 2  # capture
    self._color = "b"
    self.Capture += 1
    self.alive = 0


def inelastic_scattering(self: object) -> None:
    """
    Handles the inelastic scattering interaction of a neutron.

    Notes
    -----
    This function updates the interaction type to inelastic scattering, sets
    the neutron as alive, and increments the scattering count. It calculates
    the collision energy loss parameter and updates the neutron's energy,
    direction, and speed. The energy is reduced until thermalization.

    The function sets the color to yellow to represent inelastic scattering events.
    """

    self.intertype = 3  # scatter
    self._color = "y"
    self.Scattering += 1
    self.alive = 1
    # =====================
    # Neutron Slowing Down
    # =====================
    self.ksi = 1 + np.log((self.A - 1) / (self.A + 1)) * (self.A - 1) ** 2 / (
        2 * self.A
    )  # Collision Energy Loss Parameter
    new_E = max(
        THERMAL_E, self.energy_i * np.exp(-self.ksi)
    )  # new energy until thermalization
    self.energy_i = new_E
    self.theta_i += (
        np.pi * (1 - np.exp(-self.ksi)) * np.random.choice([-1, 1])
    )  # new direction
    self.speed_i = 1.383e7 * np.sqrt(self.energy_i)  # current neutron speed [m/s]


def pid_control(self: object) -> None:
    """
    PID control algorithm to adjust the absorption coefficient.

    Parameters
    ----------
    None

    Notes
    -----
    The function calculates the error between the current power output and the target power.
    The error is then used to calculate the proportional and integral terms of the PID
    controller. The absorption coefficient is updated by multiplying it with the proportional
    term and adding the integral term. The updated absorption coefficient is then clipped to
    the range [0, 1].

    The PID controller coefficients are set to the constants PROP_K and INT_K.
    """

    error = self.power_output / TARGET_PWR - 1
    p = 1 + PROP_K * error
    i = INT_K * error
    self.abs_coeff *= p
    self.abs_coeff += i

    self.abs_coeff = np.clip(self.abs_coeff, 0, 1)
