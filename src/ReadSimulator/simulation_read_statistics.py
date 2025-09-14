import sys
import matplotlib.pyplot as plt
import textwrap
from matplotlib.ticker import FuncFormatter


def main():
    # Ensure 5 numbers are provided as arguments
    if len(sys.argv) != 7:
        print("Usage: python bar_plot.py <num1> <num2> <num3> <num4> <num5> <num6>")
        sys.exit(1)

    try:
        # Convert arguments to integers
        values = [int(arg) for arg in sys.argv[1:]]
    except ValueError:
        print("Error: All arguments must be integers.")
        sys.exit(1)

    labels = [
        "all reads",
        "non-split reads",
        "non-split reads, no mismatches",
        "split reads",
        "split reads, no mismatches",
        "split reads, no mismatches, regions at least 5bp"
    ]
    wrapped_labels = [textwrap.fill(label, 10) for label in labels]

    # Function to format large numbers with spaces
    def format_with_spaces(number, _):
        return f'{int(number):,}'.replace(",", " ")  # Replace commas with spaces

    # Create a bar plot
    plt.figure(figsize=(12, 6))  # Adjust size for better readability
    bars = plt.bar(wrapped_labels, values, color='skyblue', edgecolor='black')

    # Add title and labels
    plt.title("Comparison of Read Types", fontsize=14)
    plt.ylabel("Frequency", fontsize=12)

    # Ensure y-axis shows full numbers (disable scientific notation)
    plt.gca().yaxis.set_major_formatter(FuncFormatter(format_with_spaces))

    # Adjust y-axis limits dynamically
    plt.ylim(0, max(values) + max(values) * 0.1)  # Add 10% padding above tallest bar

    # Add formatted value labels above each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,  # Center of the bar
            height + max(values) * 0.02,  # Slightly above the bar
            format_with_spaces(height, None),  # Format with spaces
            ha='center',  # Center horizontally
            va='bottom',  # Align text at the bottom
            fontsize=10
        )

    plt.tight_layout()
    plt.savefig("src/ReadSimulator/plots/bar_plot_read_statistics.png")


if __name__ == "__main__":
    main()