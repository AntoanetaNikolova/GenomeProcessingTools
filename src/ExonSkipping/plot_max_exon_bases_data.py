import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import textwrap

# ===== Constants =====
EXPERIMENTS = ["Homo_sapiens.GRCh37.67", "Homo_sapiens.GRCh37.75"]
DATA_PATHS = {
    "exons": "src/ExonSkipping/data/exons.csv",
    "bases": "src/ExonSkipping/data/bases.csv"
}
PLOT_PATHS = {
    "cumulative_exons": "src/ExonSkipping/plots/cumulative_plot_exons.png",
    "cumulative_bases": "src/ExonSkipping/plots/cumulative_plot_bases.png",
    "stats_exons": "src/ExonSkipping/plots/statistics_summary_by_exon.png",
    "stats_bases": "src/ExonSkipping/plots/statistics_summary_by_bases.png",
    "box_exons": "src/ExonSkipping/plots/box_plot_exons.png",
    "box_bases": "src/ExonSkipping/plots/box_plot_bases.png"
}
FIGSIZE = (10, 6)

# ===== Functions =====
def load_csv(path):
    with open(path, "r") as file:
        lines = file.read().strip().split("\n")
    return [list(map(float, line.split(","))) for line in lines]

def compute_statistics(data):
    stats = []
    for d in data:
        sorted_data = np.sort(d)
        median = np.median(sorted_data)
        lq = np.percentile(sorted_data, 25)
        uq = np.percentile(sorted_data, 75)
        iqr = uq - lq
        stats.append([median, lq, uq, iqr])
    return stats

def plot_cumulative(data, xlabel, ylabel, title, save_path):
    plt.figure(figsize=FIGSIZE)
    for i, d in enumerate(data):
        sorted_data = np.sort(d)
        cumulative_count = np.arange(1, len(sorted_data) + 1)
        plt.step(sorted_data, cumulative_count, where='post', label=f'GTF File {EXPERIMENTS[i]}')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(save_path, dpi=300)
    plt.close()

def plot_statistics_table(stat_list, title, save_path):
    df = pd.DataFrame(stat_list, columns=["Median", "LQ", "UQ", "IQR"], index=EXPERIMENTS)
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.axis('off')
    table = ax.table(cellText=df.values, rowLabels=df.index, colLabels=df.columns, cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)
    plt.title(title, fontsize=14)
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()

def plot_box_violin(data, title, save_path, y_ticks=None):
    wrapped_labels = [textwrap.fill(label, width=13) for label in EXPERIMENTS]
    log_data = [np.log10(np.array(d) + 1) for d in data]

    plt.figure(figsize=(15, 6))
    sns.violinplot(data=log_data, cut=0, inner="quartile")
    plt.xticks(range(len(wrapped_labels)), wrapped_labels, rotation=30)
    plt.title(title)

    if y_ticks is not None:
        plt.yticks(y_ticks)
        # Use the first and last tick values as y-axis limits
        plt.ylim(y_ticks[0], y_ticks[-1])

    plt.xlabel("GTF files")
    plt.ylabel("log10(Values + 1)")
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()

# ===== Main Workflow =====
exons_data = load_csv(DATA_PATHS["exons"])
bases_data = load_csv(DATA_PATHS["bases"])

# Cumulative plots
plot_cumulative(
    exons_data,
    xlabel="Skipped Exons",
    ylabel="Cumulative Count",
    title="Cumulative Distribution of Skipped Exons per ES-SE Event",
    save_path=PLOT_PATHS["cumulative_exons"]
)

plot_cumulative(
    bases_data,
    xlabel="Skipped Bases",
    ylabel="Cumulative Count",
    title="Cumulative Distribution of Skipped Bases per ES-SE Event",
    save_path=PLOT_PATHS["cumulative_bases"]
)

# Statistics tables
plot_statistics_table(
    compute_statistics(exons_data),
    title="Statistics Summary for max skipping exon",
    save_path=PLOT_PATHS["stats_exons"]
)

plot_statistics_table(
    compute_statistics(bases_data),
    title="Statistics Summary for max skipping bases",
    save_path=PLOT_PATHS["stats_bases"]
)

# Box/Violin plots (preserved)
plot_box_violin(
    exons_data,
    title="Box Plot of max skipped Exons per GTF file",
    save_path=PLOT_PATHS["box_exons"],
    y_ticks=np.arange(0.25, 1.26, step=0.1)
)

plot_box_violin(
    bases_data,
    title="Box Plot of max skipped Bases per GTF file",
    save_path=PLOT_PATHS["box_bases"]
)
