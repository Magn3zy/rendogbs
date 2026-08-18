import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FixedLocator

path_to_file = "/path/rendogbs/performance/ddr4.txt"
path_to_file2 = "/path/rendogbs/performance/lpddr5.txt"

df = pd.read_csv(path_to_file, sep=r"\s+")
df2 = pd.read_csv(path_to_file2, sep=r"\s+")

colors_ddr4 = {
    1: "#A8E6E1", 
    2:"#4DB6AC", 
    3: "#00796B"
}

colors_lpddr5 = {
    1: "#FDD9D2",
    2: "#F28E7F",
    3: "#D95F4A"
}

#print(df)
#print(df2)

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={"width_ratios": [1, 1]}, constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.15, wspace=0.05, hspace=0.05)

ax.scatter(df["#SIZE"], df["TIME"], color=colors_ddr4[3], label="DDR4", s = 25)
ax.scatter(df2["#SIZE"], df2["TIME"], color=colors_lpddr5[3], label="LPDDR5", s = 25)
ax.set_xlabel("Genome Size")

def genome_size(x, pos):
    if x == 0:
        return "0"
    elif x < 1e9:
        return f"{x/1e6:.0f} Mb"
    else:
        return f"{x/1e9:.2g} Gb"

    
#dynamicky grid na datech
ax.xaxis.set_major_locator(FixedLocator(sorted(df["#SIZE"])))
ax.xaxis.set_major_formatter(FuncFormatter(genome_size))
ax.tick_params(axis="x", labelrotation=45, labelsize=6)

#ax.xaxis.set_major_formatter(
#    FuncFormatter(lambda x, pos: "0" if x == 0 else rf"${x:.0e}".replace("e+", r"\times10^{") + "}$")
#)

#ax.set_xscale("log", base=2)
ax.grid(True, color="gray", alpha=0.25, linestyle="--")
ax.set_ylabel("Time (s)")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

ax.errorbar(df["#SIZE"], df["TIME"], yerr=df["TIMESTDEV"],
            fmt="none", ecolor=colors_ddr4[1], elinewidth=1, capsize=3, zorder=1)

ax.errorbar(df2["#SIZE"], df2["TIME"], yerr=df2["TIMESTDEV"],
            fmt="none", ecolor=colors_lpddr5[1], elinewidth=1, capsize=3, zorder=1)

ax.legend()

ax2.scatter(df["#SIZE"], df["MEMORY"], color=colors_ddr4[3], label="DDR4", s = 25)
ax2.scatter(df2["#SIZE"], df2["MEMORY"], color=colors_lpddr5[3], label="LPDDR5", s = 25)
ax2.set_xlabel("Genome Size")


#dynamicky grid na datech
ax2.xaxis.set_major_locator(FixedLocator(sorted(df["#SIZE"])))
ax2.xaxis.set_major_formatter(FuncFormatter(genome_size))
ax2.tick_params(axis="x", labelrotation=45, labelsize=6)

#ax2.set_xscale("log", base=2)
ax2.grid(True, color="gray", alpha=0.25, linestyle="--")
plt.ylabel("Memory (GB)")
ax2.set_xlim(left=0)
ax2.set_ylim(bottom=0)

ax2.errorbar(df["#SIZE"], df["MEMORY"], yerr=df["MEMSTDEV"],
            fmt="none", ecolor=colors_ddr4[1], elinewidth=1, capsize=3, zorder=1)

ax2.errorbar(df2["#SIZE"], df2["MEMORY"], yerr=df2["MEMSTDEV"],     
            fmt="none", ecolor=colors_lpddr5[1], elinewidth=1, capsize=3, zorder=1)
ax2.legend()

fig.savefig("/path/graph_dynamic.png", dpi=300)