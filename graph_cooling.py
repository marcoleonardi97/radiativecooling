import glob
import matplotlib.pyplot as plt

ids = ['T', 'Ex H', 'Ex He+', 'Ion H', 'Ion He', 'Ion He+', 'Recomb H+',
        'Recomb He+', 'Recomb He++', 'Dielect He+', 'Brehm', 'Lambda']
    
files = glob.glob("*.txt")
fls = [i for i in files if i[0] == "c"]
for i, f in enumerate(fls):
    with open(f) as data:
        rows = data.readlines()
        rows = [row.rstrip().split() for row in rows]
    T = [float(row[0]) for row in rows]
    exH = [float(row[1]) for row in rows]
    exHep = [float(row[2]) for row in rows]
    ionH= [float(row[3]) for row in rows]
    ionHe = [float(row[4]) for row in rows]
    ionHep = [float(row[5]) for row in rows]
    recombHp = [float(row[6]) for row in rows]
    recombHep = [float(row[7]) for row in rows]
    recombHepp = [float(row[8]) for row in rows]
    dielHep = [float(row[9]) for row in rows]
    brehm = [float(row[10]) for row in rows]
    lambdaf = [float(row[11]) for row in rows]

    fig = plt.figure()
    fig.set_size_inches(6,6)
    plt.title(f"Cooling function of HHe for X = {f[9:-12]}, nH = {f[16:-4]}")
    plt.ylim(-26, -20)
    plt.xlabel("Log(T) - K")
    plt.ylabel("Radiative Cooling - erg cm^3 / s")
    plt.plot(T, exH, label=ids[1], alpha=0.3)
    plt.plot(T, exHep, label=ids[2], alpha=0.3)
    plt.plot(T, ionH, label=ids[3], alpha=0.3)
    plt.plot(T, ionHe, label=ids[4], alpha=0.3)
    plt.plot(T, ionHep, label=ids[5], alpha=0.3)
    plt.plot(T, recombHp, label=ids[6], alpha=0.3)
    plt.plot(T, recombHep, label=ids[7], alpha=0.3)
    plt.plot(T, recombHepp, label=ids[8], alpha=0.3)
    plt.plot(T, dielHep, label=ids[9], alpha=0.3)
    plt.plot(T, brehm, label=ids[10], alpha=0.3)
    plt.plot(T, lambdaf, label=ids[11])
    plt.scatter(T, lambdaf, label=ids[11], alpha=0.3)
    plt.legend(loc='upper right')
    plt.savefig(f"Cooling function for X = {f[9:-12]}, nH = {f[16:-4]}.png")
    plt.show()
