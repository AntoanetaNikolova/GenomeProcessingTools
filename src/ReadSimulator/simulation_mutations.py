import matplotlib.pyplot as plt


def main():
    with open("src/ReadSimulator/output/mutations.csv", "r") as file:
        data = list(map(int, file.read().strip().split(",")))


    plt.figure(figsize=(12, 6))
    plt.hist(data, bins= len(data), edgecolor='black', alpha=0.7)
    plt.title('Number of mutations per read pair')
    plt.xlabel('Number of mutations')
    plt.ylabel('Frequency')

    plt.savefig('src/ReadSimulator/plots/mutation_distribution_plot_new.png')


if __name__ == "__main__":
    main()