import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr
from matplotlib import pyplot as plt
import argparse


def plot_correlation(experimental_data_path: str, boltz2_data_path: str):
    """
    Plots the correlation between Boltz affinity predictions and experimental IC50 values.

    Args:
        experimental_data_path (str): Path to the CSV file containing experimental IC50 data.
        boltz2_data_path (str): Path to the CSV file containing Boltz data.

    Returns:
        None
    """

    combined = pd.read_csv(experimental_data_path)
    
    combined["log_ic50"] = np.log(combined["ic50"])
    boltz2_data = pd.read_csv(
        boltz2_data_path
    )
    df = pd.merge(boltz2_data, combined, on="zinc_id", how="inner")
    df = df.dropna(subset=["ic50", "log_ic50"])

    x = df["affinity_probability_binary"]
    y = df["log_ic50"]
    spearman_corr, _ = spearmanr(x, y)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter plot
    ax.scatter(x, y, alpha=0.5)

    # Add trend line
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line = slope * x + intercept
    ax.plot(x, line, color="red", linestyle="--", alpha=0.8)

    # Customize plot
    ax.set_xlabel("Boltz Affinity (log IC50)")
    ax.set_ylabel("log IC50")
    ax.set_title(f"Boltz Affinity Binary, Spearman Correlation = {spearman_corr:.2f}")

    # Add grid
    ax.grid(True, alpha=0.3)

    # Adjust layout
    plt.tight_layout()
    plt.savefig(
        "affinity_probability_binary.png",
        dpi=600,
    )

    x = df["affinity_pred_value"]
    y = df["log_ic50"]
    spearman_corr, _ = spearmanr(x, y)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter plot
    ax.scatter(x, y, alpha=0.5)

    # Add trend line
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line = slope * x + intercept
    ax.plot(x, line, color="red", linestyle="--", alpha=0.8)

    # Customize plot
    ax.set_xlabel("Boltz Affinity (log IC50)")
    ax.set_ylabel("log IC50")
    ax.set_title(f"Boltz Affinity, Spearman Correlation = {spearman_corr:.2f}")

    # Add grid
    ax.grid(True, alpha=0.3)

    # Adjust layout
    plt.tight_layout()
    plt.savefig(
        "affinity_probability.png",
        dpi=600,
    )

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the correlation between Boltz affinity predictions and experimental IC50 values.")
    parser.add_argument("--experimental_data_path", type=str, required=True)
    parser.add_argument("--boltz2_data_path", type=str, required=True)
    args = parser.parse_args()
    plot_correlation(args.experimental_data_path, args.boltz2_data_path)