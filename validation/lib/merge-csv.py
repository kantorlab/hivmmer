import csv
import sys

print("Dataset", "Gene", "DBLength", "Evalue", "Score", "Length", "Read", sep=",")

dblength = {
    "env": 854*5226,
    "gag": 429*5183,
    "nef": 207*5188,
    "pol": 983*2956,
    "tat": 71*3361,
    "vif": 171*3937,
    "vpr": 89*3931,
    "vpu": 53*4573
}

for csv_file in sys.argv[1:]:
    dataset = csv_file.split("/")[-2]
    with open(csv_file) as f:
        for row in csv.DictReader(f):
            print(dataset,
                  row["query_name"],
                  dblength[row["query_name"]],
                  row["evalue"],
                  row["score"],
                  row["target_length"],
                  row["target_name"].partition("-")[0],
                  sep=",")

