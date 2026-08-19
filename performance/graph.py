import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FixedLocator
import numpy as np
from scipy.optimize import curve_fit


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

def func(x, a, b):
    return a * x + b

def r_squared(y_actual, y_pred):
    ss_res = np.sum((y_actual - y_pred) ** 2)
    ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
    return 1 - ss_res / ss_tot

x = np.array(df["#SIZE"])
y = np.array(df["TIME"])
x2 = np.array(df2["#SIZE"])
y2 = np.array(df2["TIME"])

popt, _ = curve_fit(func, x, y, sigma=df["TIMESTDEV"], absolute_sigma=True)
popt2, _ = curve_fit(func, x2, y2, sigma=df2["TIMESTDEV"], absolute_sigma=True)

#ax.plot(x, func(x, *popt), color=colors_ddr4[1], alpha=0.5)
#ax.plot(x2, func(x2, *popt2), color=colors_lpddr5[1], alpha=0.5)

a, b = popt
r2 = r_squared(y, func(x, *popt))
r2_2 = r_squared(y2, func(x2, *popt2))
a2, b2 = popt2

x_smooth = np.linspace(0, x.max(), 100)
x2_smooth = np.linspace(0, x2.max(), 100)

#ax.plot(x_smooth, func(x_smooth, *popt), color=colors_ddr4[1], alpha=0.5,
 #       label=f"DDR4: y = {a:.2e}x + {b:.2f}, $R^2$={r2:.3f}")
#ax.plot(x2_smooth, func(x2_smooth, *popt2), color=colors_lpddr5[1], alpha=0.5,
  #      label=f"LPDDR5: y = {a2:.2e}x + {b2:.2f}, $R^2$={r2_2:.3f}")

ax.plot(x_smooth, func(x_smooth, *popt), color=colors_ddr4[1], alpha=0.5,
        label=rf"DDR4: $y = {a/10**int(np.floor(np.log10(abs(a)))):.2f}\times10^{{{int(np.floor(np.log10(abs(a))))}}}x + {b:.2f}$, $R^2={r2:.3f}$")
ax.plot(x2_smooth, func(x2_smooth, *popt2), color=colors_lpddr5[1], alpha=0.5,
        label=rf"LPDDR5: $y = {a2/10**int(np.floor(np.log10(abs(a2)))):.2f}\times10^{{{int(np.floor(np.log10(abs(a2))))}}}x + {b2:.2f}$, $R^2={r2_2:.3f}$")

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

x3 = np.array(df["#SIZE"])
y3 = np.array(df["MEMORY"])
x4 = np.array(df2["#SIZE"])
y4 = np.array(df2["MEMORY"])

popt3, _ = curve_fit(func, x3, y3, sigma=df["MEMSTDEV"], absolute_sigma=True)  
popt4, _ = curve_fit(func, x4, y4, sigma=df2["MEMSTDEV"], absolute_sigma=True) 

a3, b3 = popt3
r2_3 = r_squared(y3, func(x3, *popt3))
r2_4 = r_squared(y4, func(x4, *popt4))
a4, b4 = popt4

x3_smooth = np.linspace(0, x3.max(), 100)
x4_smooth = np.linspace(0, x4.max(), 100)

#ax2.plot(x3_smooth, func(x3_smooth, *popt3), color=colors_ddr4[1], linewidth=1.5, alpha=0.4,
        #label=f"DDR4: y = {a3:.2e}x + {b3:.2f}, $R^2$={r2_3:.3f}")
#ax2.plot(x4_smooth, func(x4_smooth, *popt4), color=colors_lpddr5[1], alpha=0.4, linewidth=1.5,
        #label=f"LPDDR5: y = {a4:.2e}x + {b4:.2f}, $R^2$={r2_4:.3f}")

ax2.plot(x3_smooth, func(x3_smooth, *popt3), color=colors_ddr4[1], linewidth=1.5, alpha=0.4,
        label=rf"DDR4: $y = {a3/10**int(np.floor(np.log10(abs(a3)))):.2f}\times10^{{{int(np.floor(np.log10(abs(a3))))}}}x + {b3:.2f}$, $R^2={r2_3:.3f}$")
ax2.plot(x4_smooth, func(x4_smooth, *popt4), color=colors_lpddr5[1], alpha=0.4, linewidth=1.5,
        label=rf"LPDDR5: $y = {a4/10**int(np.floor(np.log10(abs(a4)))):.2f}\times10^{{{int(np.floor(np.log10(abs(a4))))}}}x + {b4:.2f}$, $R^2={r2_4:.3f}$")


ax.errorbar(df["#SIZE"], df["TIME"], yerr=df["TIMESTDEV"],
            fmt="none", ecolor=colors_ddr4[1], elinewidth=1, capsize=3, zorder=1)

ax.errorbar(df2["#SIZE"], df2["TIME"], yerr=df2["TIMESTDEV"],
            fmt="none", ecolor=colors_lpddr5[1], elinewidth=1, capsize=3, zorder=1)


ax2.scatter(df["#SIZE"], df["MEMORY"], color=colors_ddr4[3], label="DDR4", s = 25)
ax2.scatter(df2["#SIZE"], df2["MEMORY"], color=colors_lpddr5[3], label="LPDDR5", s = 25)
ax2.set_xlabel("Genome Size")

ax.legend(fontsize=8)

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
ax2.legend(fontsize=8)

#plt.show()
fig.savefig("/path/rendogbs/performance/graph_function.png", dpi=300)