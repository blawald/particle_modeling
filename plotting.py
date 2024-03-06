import os
import matplotlib.pyplot as plt
directories = [f for f in os.listdir('.') if (os.path.isdir(f) and f!=".vim" and f!=".git")]
for current_dir in directories:
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax_extra = fig.add_subplot(212)
    with open(current_dir + "//Energy_data.xyz") as f:
        lines = f.readlines()
    names = lines[0].split('\t')[:-1]
    data = [[] for i in range(len(names))]
    for i in range(1, len(lines)-1):
        for j in range(len(names)):
            temp = lines[i].split('\t')[:-1]
            data[j].append(temp[j])
    
    ax.plot(list(map(float, data[0])), list(map(float, data[1])), label=names[1])
    ax.plot(list(map(float, data[0])), list(map(float, data[2])), label=names[2])
    ax.plot(list(map(float, data[0])), list(map(float, data[3])), label=names[3])
    ax.grid()
    ax.set_ylabel("Energy")
    ax.legend()
    ax_extra.plot(list(map(float, data[0])), list(map(float, data[1])), label = names[1])
    ax_extra.legend()
    ax_extra.grid()
    ax_extra.set_xlabel("Number of iterations")
    ax_extra.set_ylabel("Energy")
    plt.savefig(current_dir + "//Enery_analysis.png")
    with open(current_dir + "//v_data.xyz") as f:
        lines = f.readlines()
    fig = plt.figure()
    ax = fig.add_subplot()
    t = float(int(lines[0])+2)
    print(t)
    i = int(len(lines)/(t)) - 1
    print(i)
    data = []
    print(current_dir)
    for j in range(int(t*i+2), int(t*i+t)):
        data.append(float(lines[j].split('\t')[1]))
    plt.clf()
    plt.hist(data, bins=35)
    plt.savefig(current_dir+"//Velocity_x_hist")
