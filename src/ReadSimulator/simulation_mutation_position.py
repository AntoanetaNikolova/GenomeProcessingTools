import matplotlib.pyplot as plt


def main():
    with open("src/ReadSimulator/output/mutationPosition.csv", "r") as file:
        data = list(map(int, file.read().strip().split(",")))


    # Histogram
    plt.hist(data, bins= len(data), edgecolor='black', alpha=0.7)
    plt.title('Fragment Lengths Distribution')
    plt.xlabel('Fragment lengths')
    plt.ylabel('Frequency')

    plt.savefig('src/ReadSimulator/plots/mutationPosition.png')

if __name__ == "__main__":
    main()