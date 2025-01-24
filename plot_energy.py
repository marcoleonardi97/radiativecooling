import glob
import matplotlib.pyplot as plt

files = glob.glob("*.txt")
fls = [i for i in files if i[0] == "e"]
for f in fls:
    with open(f) as hdul:
        rows = hdul.readlines()
        rows = [row.rstrip().split() for row in rows]
    time = [float(row[0]) for row in rows]
    temperature = [float(row[1]) for row in rows]
    fig = plt.figure()
    fig.set_size_inches(5,5)
    plt.title(f"Time evolution of temperature at nH: {f[23:-4]}")
    plt.yscale("log")
    plt.xscale("log")
    plt.scatter(time, temperature, c="red")
    plt.plot(time, temperature, alpha=0.5, c="red")
    plt.xlabel("Log of Cooling Time")
    plt.ylabel("Log T - Kelvin")
    plt.savefig(f"temp_evo_nH_{f[23:-4]}.png")
    plt.show()
