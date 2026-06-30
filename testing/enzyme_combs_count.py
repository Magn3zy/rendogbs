import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

input_file = "/Users/eliskakorbova/Desktop/biopython/rendogbs/testing/enzyme_combs_counts.csv"

# Load data
df = pd.read_csv(input_file)

# Sort by count
df = df.sort_values("count")

# Color palettes
palette_no = {
    1: "#A8E6E1",   # light green-turquoise
    2: "#4DB6AC",   # medium turquoise
    3: "#00796B"    # dark teal
}

palette_yes = {
    1: "#FDD9D2",
    2: "#F28E7F",
    3: "#D95F4A"
}

# Assign colors according to count AND methylation sensitivity
colors = []

for _, row in df.iterrows():
    count = row["count"]
    methyl = row["methylation_sensitive"].strip().lower()

    if methyl == "yes":
        colors.append(palette_yes[count])
    else:
        colors.append(palette_no[count])

# Plot
fig, ax = plt.subplots(figsize=(10, 8))

bars = ax.barh(
    df["combination"],
    df["count"],
    color=colors
)

# Values at the end of bars
for bar in bars:
    width = bar.get_width()
    ax.text(
        width + 0.05,
        bar.get_y() + bar.get_height() / 2,
        str(int(width)),
        va="center",
        fontsize=10
    )

# Axes
ax.set_xlabel("Counts")
ax.set_ylabel("Restriction enzyme combination")
ax.set_title("Restriction enzyme combinations")

ax.xaxis.set_major_locator(MultipleLocator(1))
ax.set_ylim(-0.5, len(df) - 0.5)

# Clean style
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", alpha=0.25)
ax.set_axisbelow(True)

# Legend
from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor="#4DB6AC", label="No methylation-sensitive enzyme"),
    Patch(facecolor="#F28E7F", label="Contains methylation-sensitive enzyme")
]

ax.legend(
    handles=legend_elements,
    frameon=False,
    loc="lower right"
)

plt.tight_layout()

# Save
output = input_file.replace(".csv", "_barplot.png")
plt.savefig(output, dpi=300, bbox_inches="tight")

plt.show()