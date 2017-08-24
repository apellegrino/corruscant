import matplotlib.pyplot as plt
import numpy as np
import richardsplot

def line_parse(line, delimiter=','):
    return filter(None, line.strip('\r\n').split(delimiter)) if line else None

plt.figure(figsize=(8,8))

n_arr, avg_arr, err_arr = [], [], []
with open("results.csv", 'r') as f:
    header = f.readline()

    for line in f:
        n, avg, err = line_parse(line)
        n = int(n)
        avg, err = float(avg), float(err)

        n_arr.append(n)
        avg_arr.append(avg)
        err_arr.append(err)

#plt.errorbar(n_arr, avg_arr, yerr=err_arr, capsize=2.)
plt.plot(n_arr, avg_arr, 'o', color='blue', ls='-')
plt.ylim((0,9.0))
plt.xlabel("Number of fields")
plt.ylabel("Average computation time (sec), n = 5")
plt.tight_layout()
plt.savefig("field_times.png", dpi=200)
