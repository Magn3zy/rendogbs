import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator


input_file = "/Users/eliskakorbova/Desktop/biopython/rendogbs/testing/enzyme_combs_counts.csv"

df = pd.read_csv(input_file)

df = df.sort_values("count")

norm = Normalize(df["count"].min(), df["count"].max())
palette = {
    1: "#A8E6E1",
    2: "#4DB6AC",
    3: "#00796B"
}

colors = [palette[c] for c in df["count"]]

fig, ax = plt.subplots(figsize=(10,8))

bars = ax.barh(
    df["combination"],
    df["count"],
    color=colors
)

for bar in bars:
    width = bar.get_width()
    ax.text(width + 0.05,
            bar.get_y() + bar.get_height()/2,
            str(int(width)),
            va='center')

ax.set_xlabel("Counts")
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.set_ylabel("Restriction enzyme combination")
ax.set_ylim(-0.5, len(df) - 0.5)
ax.set_title("RE combinations counts")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()


output = input_file.replace(".csv", "_barplot.png")
plt.savefig(output, dpi=300)

plt.show()