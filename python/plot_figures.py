#!/usr/bin/python

import os

input = "../out"
outdir = "output"
os.system("mkdir -p %s" % outdir)

outprop = "%s/prop" % outdir
outmass = "%s/mass" % outdir
title = "$C_8$"

cmds = []
cmds.append( "../../nucnet-tools-code/examples/analysis/print_properties" )
cmds.append( "../../nucnet-tools-code/examples/analysis/print_mass_fractions_in_zones" )

os.system("python make_analysis.py")
os.system("%s %s time t9 rho > %s" % (cmds[0], input, outprop))
os.system("%s %s n h1g h2 he2c li3c be4c b5c c6c n7c o8c o8r > %s " \
          % (cmds[1], input, outmass))

names = ['$O$', '$C$', '$10^3CO$', 
         '$10^{10}C_2$', '$10^{16}C_3$', '$10^{15}C_4$',
         '$10^{17}C_5$', '$10^{15}C_6$',
         '$10^{18}C_7$', '$10^{23}C_8^C$', '$10^{16}C_8^R$']

import matplotlib.pyplot as plt

time = []
t9 = []
rho = []

for line in open("%s" % outprop, "r").readlines():

    line = line[:-1].split()

    time.append(float(line[1]))
    t9.append(line[2])
    rho.append(line[3])

time0 = time[0]

for i in range(len(time)):
    time[i] -= time0 

masses = []

for line in open("%s" % outmass, "r").readlines():

    line = line[:-1].split()
    masses.append(line)

xx = []

for i in range(1,len(masses[0])):
    x = []
    for j in range(len(masses)):
        x.append(float(masses[j][i]))
    xx.append(x)

factors = [2, 2, 3, 4, 5, 6, 7, 8, 8]
powers = [3, 10, 16, 15, 17, 15, 18, 23, 16]

for i in range(2,len(xx)):
     for j in range(len(xx[i])): 
         factor = 10 ** powers[i-2] / factors[i-2]
         xx[i][j] *= factor

for i in range(len(xx)-5):
    plt.loglog(time,xx[i],label=names[i])

for i in range(len(xx)-5, len(xx)):
    plt.loglog(time,xx[i],'--', label=names[i])

plt.axis([1e4,2e8,1e-8,1e1])
plt.title(title)
plt.legend(loc = "center left")
plt.savefig("%s/mass.png" % outdir)
#plt.show()
