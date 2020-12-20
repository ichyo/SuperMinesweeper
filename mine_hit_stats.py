import csv
import sys
import numpy as np
import matplotlib.pyplot as plt

path = sys.argv[1]
values = []
with open(path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['reason'] == 'win' or row['reason'] == 'maybe_win':
            values.append(int(row['mine_hit']))

values = np.array(values)
print(values.mean())
print(values.std())
