import os
import numpy
import matplotlib.pyplot as plt
directories = [f for f in os.listdir('.') if (os.path.isdir(f) and f!=".vim" and f!=".git" and f!="env")]
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
    data = []
    lower_border = int(1*i/2.0)
    if lower_border > i: lower_border = i
    for j in range(lower_border, i):
        for k in range(2, int(t)):
            data.append(float(lines[int(j*t + k)].split('\t')[1]))
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.hist(data, bins=150)
    ax.grid()
    ax.set_xlabel("V_x")
    ax.set_ylabel("Number of particles")
    plt.savefig(current_dir+"//Velocity_x_hist_time")
    fig = plt.figure()
    ax = fig.add_subplot()
    data = list(map(lambda x: x*x, data))
    plt.hist(data, bins=150)
    ax.set_yscale('log')
    ax.grid()
    ax.set_xlabel("V_x squared")
    ax.set_ylabel("Number of particles")
    plt.savefig(current_dir+"//Velocity_x_hist_time_linear")
    f = lambda x: numpy.log(-x*x)
    data = []
    temp
    for j in range(0, i):
        temp = 0
        for k in range(2, int(t)):
            temp+=float(lines[int(j*t + k)].split('\t')[1])
        data.append(temp)
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.plot([100*j for j in range(0, i)], data)
    plt.savefig(current_dir+"//P_analysis")
    
    fig = plt.figure()
    ax = fig.add_subplot()
    data = list(map(lambda x: x*x, data))
    plt.hist(data, bins=150)
    ax.set_yscale('log')
    ax.grid()
    ax.set_xlabel("V_x squared")
    ax.set_ylabel("Number of particles")
    f = lambda x: numpy.log(-x*x)
    tmp = 0.0
    data = []
    higher_border = lower_border + 150
    tmp1 = 0
    if higher_border > i: higher_border = i
    for j in range(lower_border, higher_border):
        data1 = []
        for k in range(2, int(t)):
            tmp1 = 0
            tmp = float(lines[int(j*t + k)].split('\t')[1])
            tmp -= float(lines[int(lower_border*t + k)].split('\t')[1])
            tmp = tmp*tmp
            tmp1 +=tmp
            tmp = float(lines[int(j*t + k)].split('\t')[2])
            tmp -= float(lines[int(lower_border*t + k)].split('\t')[2])
            tmp = tmp*tmp
            tmp1 +=tmp
            tmp = float(lines[int(j*t + k)].split('\t')[3])
            tmp -= float(lines[int(lower_border*t + k)].split('\t')[3])
            tmp = tmp*tmp
            tmp1 +=tmp
            data1.append(tmp1)
        data.append(numpy.average(data1))
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.scatter([100*j for j in range(lower_border, higher_border)], data)
    plt.savefig(current_dir+"//<delta_r^2>(t)")
    plt.close()
