import toml
import matplotlib.pyplot as plt
from math import ceil

DPI = 1000  # Define the resolution for saved plots


def plot_results(n_values, variable_name, values, plot_filename, mark_min=False):
    """
    Plots a graph of the given values against the number of fins and saves the plot.

    Parameters:
    n_values (iterable): The x-axis values (number of fins).
    variable_name (str): The label for the y-axis (variable name).
    values (iterable): The y-axis values (computed values for the variable).
    plot_filename (str): The file path to save the plot.
    mark_min (bool): Whether to mark the minimum value on the plot.
    """
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(n_values, values)
    ax1.set_xlabel("Number of fins (n)")
    ax1.set_ylabel(variable_name)
    ax1.grid()

    if mark_min:
        # Mark the minimum value and annotate it
        min_value = min(values)
        min_index = values.index(min_value)
        ax1.plot(n_values[min_index], min_value, "o", color="#F97306")
        ax1.annotate(
            f"Min: {min_value:.2f}",
            (n_values[min_index], min_value),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            color="#F97306",
        )
        print(
            f"Minimum {variable_name} value: {min_value} at n = {n_values[min_index]}"
        )

    fig.savefig(plot_filename, dpi=DPI)
    plt.close(fig)


class HeatSink:
    """
    A class representing a heat sink model for calculating various heat-related parameters.
    """

    def __init__(self, config):
        """
        Initializes the HeatSink object with the parameters from the provided configuration.

        Parameters:
        config (dict): Configuration values from a TOML file.
        """
        # Load configuration values
        self.k_silicon = config["silicon"]["k"]
        self.k_water = config["water"]["k"]
        self.mu_water = config["water"]["mu"]
        self.c_p_water = config["water"]["cp"]
        self.rho_water = config["water"]["rho"]
        self.q = config["power"]
        self.w_chip = config["dimensions"]["w_chip"]
        self.l_chip = config["dimensions"]["l_chip"]
        self.h_fins = config["dimensions"]["h_fins"]
        self.t = config["dimensions"]["t"]
        self.dp_max = config["dp_max"]
        self.t_in = config["t_in"]
        self.nu = config["nu"]
        self.min_w = config["min_wc"]

    def compute_flow_rate(self, w_c, n):
        """
        Computes the flow rate based on the fin width and number of fins.

        Parameters:
        w_c (float): The width of the coolant channel.
        n (int): The number of fins.

        Returns:
        float: The computed flow rate.
        """
        return (
            self.dp_max
            * (1 - 0.63 * w_c / self.h_fins)
            * (n + 1)
            * (w_c**3)
            * self.h_fins
            / (12 * self.mu_water * self.l_chip)
        )

    def compute_chip_conductance(self, w_c=None, n=None):
        """
        Computes the conductance of the chip based on its dimensions and thermal properties.

        Returns:
        float: The chip conductance.
        """
        a_chip = self.w_chip * self.l_chip
        return self.t / (self.k_silicon * a_chip)

    def compute_convective_resistance(self, w_c, n=None):
        """
        Computes the convective resistance between the chip and coolant.

        Parameters:
        w_c (float): The width of the coolant channel.
        n (int): The number of fins (optional).

        Returns:
        float: The convective thermal resistance.
        """
        a_w = (w_c + 2 * self.h_fins) * self.w_chip * self.l_chip / (2 * w_c)
        htc = self.k_water * self.nu / w_c
        return 1 / (htc * a_w)

    def compute_sensible_heat_resistance(self, w_c, n):
        """
        Computes the sensible heat resistance based on the flow rate.

        Parameters:
        w_c (float): The width of the coolant channel.
        n (int): The number of fins.

        Returns:
        float: The sensible heat resistance.
        """
        f = self.compute_flow_rate(w_c, n)
        return 1 / (self.rho_water * self.c_p_water * f)

    def compute_total_resistance(self, w_c, n):
        """
        Computes the total thermal resistance considering conductive, convective, and heat resistances.

        Parameters:
        w_c (float): The width of the coolant channel.
        n (int): The number of fins.

        Returns:
        float: The total thermal resistance.
        """
        r_heat = self.compute_sensible_heat_resistance(w_c, n)
        r_cond = self.compute_chip_conductance()
        r_conv = self.compute_convective_resistance(w_c, n)
        return r_cond + r_conv + r_heat

    def compute_surface_temperature(self, w_c, n):
        """
        Computes the surface temperature of the chip based on the total thermal resistance.

        Parameters:
        w_c (float): The width of the coolant channel.
        n (int): The number of fins.

        Returns:
        float: The surface temperature of the chip.
        """
        r_total = self.compute_total_resistance(w_c, n)
        return self.t_in + r_total * self.q * self.w_chip * self.l_chip

    def compute_wc(self, n):
        """
        Computes the coolant channel width based on the number of fins.

        Parameters:
        n (int): The number of fins.

        Returns:
        float: The coolant channel width.
        """
        return self.w_chip / (2 * n + 1)

    def compute_all_results(self, n_min, n_max, method):
        """
        Computes the results for a range of fins based on the provided method.

        Parameters:
        n_min (int): The minimum number of fins.
        n_max (int): The maximum number of fins.
        method (function): The method to compute the results (e.g., heat resistance).

        Returns:
        tuple: A tuple of n_values (x-values) and results (y-values).
        """
        n_values = range(n_min, n_max)
        results = []

        for n in n_values:
            w_c = self.compute_wc(n)
            result = method(w_c, n)
            results.append(result)

        return n_values, results


# Load configuration from TOML file
config = toml.load("config.toml")

# Create the HeatSink object with the loaded configuration
heat_sink = HeatSink(config)

# Exploration range for fins
n_max = ceil((heat_sink.w_chip / heat_sink.min_w - 1) / 2)
n_min = ceil((heat_sink.w_chip / heat_sink.h_fins - 1) * 0.5)

# Compute and plot various results
n_values, r_heat_values = heat_sink.compute_all_results(
    n_min, n_max, heat_sink.compute_sensible_heat_resistance
)
plot_results(
    n_values, r"$R_{heat} \left[\frac{K}{W}\right]$", r_heat_values, "images/r_heat.png"
)

n_values, r_conv_values = heat_sink.compute_all_results(
    n_min, n_max, heat_sink.compute_convective_resistance
)
plot_results(
    n_values, r"$R_{conv} \left[\frac{K}{W}\right]$", r_conv_values, "images/r_conv.png"
)

n_values, r_tot_values = heat_sink.compute_all_results(
    n_min, n_max, heat_sink.compute_total_resistance
)
plot_results(
    n_values,
    r"$R_{tot} \left[\frac{K}{W}\right]$",
    r_tot_values,
    "images/r_tot.png",
    True,
)

n_values, surf_temp_values = heat_sink.compute_all_results(
    n_min, n_max, heat_sink.compute_surface_temperature
)
plot_results(
    n_values,
    r"$T_{chip} \left[K\right]$",
    surf_temp_values,
    "images/temp_chip.png",
    True,
)


# Additional plots
def vary_and_plot(
    n_min, n_max, param_values, heatsink, method, ylabel, file, param_name
):
    """
    Varies a parameter (like viscosity or pressure drop) and plots the results.

    Parameters:
    n_min (int): The minimum number of fins.
    n_max (int): The maximum number of fins.
    param_values (list): List of values for the parameter to vary.
    heatsink (HeatSink): The HeatSink object.
    method (function): The method to compute the results (e.g., flow rate).
    ylabel (str): The label for the y-axis.
    file (str): The filename to save the plot.
    param_name (str): The name of the parameter being varied.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    for param in param_values:
        setattr(heatsink, param_name, param)
        n_values, results = heatsink.compute_all_results(n_min, n_max, method)
        ax.plot(
            n_values, results, label=f"{param_name}={param}"
        )  # Add a label for legend

    ax.set_xlabel("Number of fins (n)")
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.legend(title=f"{param_name} variations")
    plt.savefig(file, dpi=DPI)
    plt.close(fig)


# Plot different resistances
fig, ax1 = plt.subplots(figsize=(10, 6))
n_values, r_cond_values = heat_sink.compute_all_results(
    n_min, n_max, heat_sink.compute_chip_conductance
)
ax1.plot(n_values, r_heat_values, label="$R_{heat}$")
ax1.plot(n_values, r_conv_values, label="$R_{conv}$")
ax1.plot(n_values, r_cond_values, label="$R_{cond}$")
ax1.plot(n_values, r_tot_values, label="$R_{tot}$", linewidth=2)
ax1.set_xlabel("Number of fins (n)")
ax1.set_ylabel("Resistance [K/W]")
ax1.grid()
ax1.legend()
fig.savefig("images/resistance_comparison.png", dpi=DPI)

# Vary viscosity and plot the results
mu_values = [1e-3, 2e-3, 3e-3]
vary_and_plot(
    n_min,
    n_max,
    mu_values,
    heat_sink,
    heat_sink.compute_flow_rate,
    r"Flow rate [$\frac{m^3}{s}$]",
    "images/mu_flowrate.png",
    "mu_water",
)

# Vary dp_max and plot the results
dp_max_values = [1e4, 5e4, 1e5]
vary_and_plot(
    n_min,
    n_max,
    dp_max_values,
    heat_sink,
    heat_sink.compute_flow_rate,
    r"Flow rate [$\frac{m^3}{s}$]",
    "images/dp_flowrate.png",
    "dp_max",
)

# Plot different resistances ad different dps
fig, ax1 = plt.subplots(figsize=(10, 6))
heatsink = HeatSink(config)
n_values, r_tot_values_1 = heat_sink.compute_all_results(
    n_min, n_max, heatsink.compute_total_resistance
)
ax1.plot(
    n_values, r_tot_values_1, label=f"$\Delta P = {heatsink.dp_max} \mathrm{{Pa}}$"
)
heatsink.dp_max = 1e5
n_values, r_tot_values_2 = heat_sink.compute_all_results(
    n_min, n_max, heatsink.compute_total_resistance
)
ax1.plot(
    n_values, r_tot_values_2, label=f"$\Delta P = {heatsink.dp_max} \mathrm{{Pa}}$"
)
ax1.set_xlabel("Number of fins (n)")
ax1.set_ylabel("Resistance [K/W]")
ax1.grid()
ax1.legend()
fig.savefig("images/resistance_comparison_dp.png", dpi=DPI)

# Plot different resistances ad different mus
fig, ax1 = plt.subplots(figsize=(10, 6))
heatsink = HeatSink(config)
n_values, r_tot_values_1 = heat_sink.compute_all_results(
    n_min, n_max, heatsink.compute_total_resistance
)
ax1.plot(
    n_values,
    r_tot_values_1,
    label=f"$\mu_{{water}} = {heatsink.mu_water} \mathrm{{Pa s}}$",
)
heatsink.mu_water = 5e-4
n_values, r_tot_values_1 = heat_sink.compute_all_results(
    n_min, n_max, heatsink.compute_total_resistance
)
ax1.plot(
    n_values,
    r_tot_values_1,
    label=f"$\mu_{{water}} = {heatsink.mu_water} \mathrm{{Pa s}}$",
)
ax1.set_xlabel("Number of fins (n)")
ax1.set_ylabel("Resistance [K/W]")
ax1.grid()
ax1.legend()
fig.savefig("images/resistance_comparison_mu.png", dpi=DPI)
