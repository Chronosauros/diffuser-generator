"""Shared chart style (dataviz reference palette, light surface, large text for reading on a phone)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#d9d8d4"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
FLAT = "#8f8e89"
BAND = "#eeede9"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": GRID, "grid.color": GRID, "grid.linewidth": 0.8,
    "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "legend.fontsize": 12,
    "legend.frameon": False, "lines.linewidth": 2.0, "axes.spines.top": False, "axes.spines.right": False,
})


def polar_axes(fig, pos, rmin=-40):
    ax = fig.add_subplot(pos, projection="polar")
    ax.set_thetamin(-90); ax.set_thetamax(90)
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)       # +theta towards +x, to the right
    ax.set_rlim(rmin, 6)
    ax.set_rticks([-30, -20, -10, 0])
    ax.set_yticklabels(["-30", "-20", "-10", "0"], fontsize=9, color=INK2)
    ax.set_rlabel_position(80)
    ax.set_xticks([-1.5708, -0.7854, 0, 0.7854, 1.5708])
    ax.set_xticklabels(["-90°", "-45°", "0°", "45°", "90°"], fontsize=10)
    ax.grid(True)
    return ax
