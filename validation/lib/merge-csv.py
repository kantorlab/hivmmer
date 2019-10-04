import csv
import sys

print("Dataset", "Gene", "Evalue", sep=",")

for csv_file in sys.argv[1:]:
    dataset = csv_file.split("/")[-2]
    with open(csv_file) as f:
        for row in csv.DictReader(f):
            print(dataset, row["target_name"], row["evalue"], sep=",")

