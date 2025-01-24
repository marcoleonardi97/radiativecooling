import glob
import matplotlib.pyplot as plt

files = glob.glob("*.txt")
fls = [i for i in files if i[0] == "i"]
for i, f in enumerate(fls):
    with open(f) as data:
        rows = data.readlines()
        rows = [row.rstrip().split() for row in rows]
    T = [float(row[0]) for row in rows]
    nH0 = [float(row[1]) for row in rows]
    nHp = [float(row[2]) for row in rows]
    nHe0 = [float(row[3]) for row in rows]
    nHep = [float(row[4]) for row in rows]
    nHepp = [float(row[5]) for row in rows]
    n_e = [float(row[6]) for row in rows]
    tot_H = [i+j for i,j in zip(nH0, nHp)]
    tot_He = [i+j+k for i,j,k in zip(nHe0, nHep, nHepp)]

    if tot_He[0] == 0: # don't plot for X = 1 or div/0
        continue
    else:

        fig = plt.figure()
        fig.set_size_inches(6,6)
        plt.title(f"Ionization state for X = {f[12:-12]}, nH = {f[19:-4]}")
        plt.xscale("log")
        plt.xlabel("K log")
        plt.ylabel("cm^3/s")
        plt.plot(T, [n/t for n,t in zip(nH0,tot_H)], label="Neutral Hydrogen")
        plt.plot(T, [n/t for n,t in zip(nHp,tot_H)], label="Ionized Hydrogen")
        plt.plot(T, [n/t for n,t in zip(nHe0,tot_He)], label="Neutral Helium")
        plt.plot(T, [n/t for n,t in zip(nHep,tot_He)], c="blue", label="Once Ionized Helium")
        plt.plot(T, [n/t for n,t in zip(nHepp,tot_He)], label="Doubly Ionized Helium")
        plt.plot(T, n_e, label="Free Electrons")
        plt.legend()
        plt.savefig(f"Ionization state for X = {f[12:-12]}, nH = {f[19:-4]}.png")
        plt.show()
