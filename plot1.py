import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull

####      Author: Aldo Abarca Ortega      #####
#### Assingment 2019. PART 1 (C. thermo.) #####
#######---Simulations on materials---##########

### Indications ##
## 1.- Units are K for temperature and J/mol for Gibbs free energy
## 2.- Stable crystal lattice for Al at 298K: fcc
## 3.- Stable crystal lattice for Zn at 298K: hcp

R = 8.314 #J/(mol*K)

# Range of temperature from 298 to 1000K
ti = 298
tf = 1000 #1000
esp = (tf-ti)/6
tem = np.linspace(ti,tf,esp)

row_num = len(tem)

sep = 100

ggfc = np.zeros((sep,len(tem)))
gghc = np.zeros((sep,len(tem)))
gglq = np.zeros((sep,len(tem)))
xxx = np.zeros((sep,len(tem)))

for x in range(1,sep):

    for i in range(len(tem)):

        lfcc = []
        lhcp = []
        lliq = []

        ## Al
        if tem[i] >= 298.00 and tem[i] < 700:

            #FCC_A1
            alfcc = -7976.15 + 137.093038 * tem[i] - 24.3671976 * tem[i] * np.log(tem[i]) - 1.884662 * 10 ** (-3) * (
                tem[i]) ** 2 - 0.877664 * 10 ** (-6) * (tem[i]) ** 3 + 74092 * (tem[i]) ** (-1)
            #HCP_A3
            alhcp = -2495.15 + 135.293038 * tem[i] - 24.3671976 * tem[i] * np.log(tem[i]) - 1.884662 * 10 ** (-3) * tem[
                i] ** 2 - 0.877664 * 10 ** (-6) * tem[i] ** 3 + 74092 * tem[i] ** (-1)
            #LIQUID
            alliq = 3028.895 + 125.251188 * tem[i] - 24.3671976 * tem[i] * np.log(tem[i]) - 1.884662 * 10 ** (-3) * tem[
                i] ** 2 - 0.877664 * 10 ** (-6) * tem[i] ** 3 + 74092 * tem[i] ** (-1) + 79.34 * 10 ** (-21) * tem[
                        i] ** 7

        elif tem[i] >= 700 and tem[i] < 933.473:

            #FCC_A1
            alfcc = -11276.24 + 223.048446 * tem[i] - 38.5844296 * tem[i] * np.log(tem[i]) + 18.531982 * 10 ** (-3) * \
                    tem[i] ** 2 - 5.764227 * 10 ** (-6) * tem[i] ** 3 + 74092 * tem[i] ** (-1)
            #HCP_A3
            alhcp = -5795.24 + 221.248446 * tem[i] - 38.5844296 * tem[i] * np.log(tem[i]) + 18.531982 * 10 ** (-3) * \
                    tem[i] ** 2 - 5.764227 * 10 ** (-6) * tem[i] ** 3 + 74092 * tem[i] ** (-1)
            #LIQUID
            alliq = -271.194 + 211.206596 * tem[i] - 38.5844296 * tem[i] * np.log(tem[i]) + 18.531982 * 10 ** (-3) * \
                    tem[i] ** 2 - 5.764227 * 10 ** (-6) * tem[i] ** 3 + 74092 * tem[i] ** (-1) + 79.34 * 10 ** (-21) * \
                    tem[i] ** (7)

        elif tem[i] >= 933.473:

            #FCC_A1
            alfcc = -11278.361 + 188.684136 * tem[i] - 31.748192 * tem[i] * np.log(tem[i]) - 1230.622 * 10 ** (25) * (
            tem[i]) ** ( -9)
            #HCP_A3
            alhcp = -5797.361+186.884136*tem[i]-31.748192*tem[i]*np.log(tem[i])-1230.622*10**(25)*tem[i]**(-9)
            #LIQUID
            alliq = -795.991+177.430209*tem[i]-31.748192*tem[i]*np.log(tem[i])

        ## ZN
        if tem[i] >= 298.00 and tem[i] < 692.677:

            #HCP_ZN
            znhcp = -7285.787 + 118.470069 * tem[i] - 23.701314 * tem[i] * np.log(tem[i]) - 1.712034 * 10 ** (-3) * (
                tem[i]) ** 2 - 1.264963 * 10 ** (-6) * (tem[i]) ** 3
            #FCC_A1
            znfcc = -4315.967 + 116.900389 * tem[i] - 23.701314 * tem[i] * np.log(tem[i]) - 1.712034 * 10 ** (-3) * tem[
                i] ** 2 - 1.264963 * 10 ** (-6) * tem[i] ** 3
            #LIQUID
            znliq = -128.565 + 108.177019 * tem[i] - 23.701314 * tem[i] * np.log(tem[i]) - 1.712034 * 10 ** (-3) * tem[
                i] ** 2 - 1.264963 * 10 ** (-6) * tem[i] ** 3 - 358.949 * 10 ** (-21) * tem[i] ** 7

        elif tem[i] >= 692.677 and tem[i] < 1700:

            #HCP_ZN
            znhcp = -11070.546 + 172.345644 * tem[i] - 31.38 * tem[i] * np.log(tem[i]) + 470.47 * 10 ** 24 * tem[i] ** (
                -9)
            #FCC_A1
            znfcc = -8100.726+170.775964*tem[i]-31.38*tem[i]*np.log(tem[i])+470.47*10**24*tem[i]**(-9)
            #LIQUID
            znliq = -3620.385+161.60844*tem[i]-31.38*tem[i]*np.log(tem[i])

        lfcc.append(7297.48+0.47512*tem[i])
        lfcc.append(6612.88-4.59110*tem[i])
        lfcc.append(-3097.19+3.30635*tem[i])

        lhcp.append(18820.95-8.95255*tem[i])
        lhcp.append(1.0*10**(-6)+0.00*tem[i])
        lhcp.append(1.0*10**(-6)+0.00*tem[i])
        lhcp.append(-702.79+0.00*tem[i])

        lliq.append(10465.55-3.39259*tem[i])

        xzn = x/sep
        xal = 1 - xzn
        xdf = xal - xzn

        gfcc = xal * alfcc + xzn * znfcc + R * tem[i] * (xal * np.log(xal) + xzn * np.log(xzn)) + xal * xzn * (
                    lfcc[0] + lfcc[1] * xdf + lfcc[2] * xdf ** 2)


        ghcp = xal * alhcp + xzn * znhcp + R * tem[i] * (xal * np.log(xal) + xzn * np.log(xzn)) + xal * xzn * (
                    lhcp[0] + lhcp[1] * xdf + lhcp[2] * xdf ** 2 + lhcp[3] * xdf ** 3)

        gliq = xal*alliq+xzn*znliq + R*tem[i]*(xal*np.log(xal)+xzn*np.log(xzn))+xal*xzn*(lliq[0])

        xxx[x][i]  = xzn
        ggfc[x][i] = gfcc
        gghc[x][i] = ghcp
        gglq[x][i] = gliq


xfin = np.zeros((len(np.transpose(xxx[:,1]))+1,len(tem)+1))

for i in range(len(xxx)):
    xfin[i, 0] = xxx[i, 1]


res1 = 0
for j in range(len(tem)):
    for i in range(len(gglq)):

        if ggfc[i][j] > gglq[i][j] and gghc[i][j] > gglq[i][j]:
            res1 = gglq[i][j]

        elif gglq[i][j] > gghc[i][j] and ggfc[i][j] > gghc[i][j]:
            res1 = gghc[i][j]

        elif gglq[i][j] > ggfc[i][j] and gghc[i][j] > ggfc[i][j]:
            res1 = ggfc[i][j]

        xfin[i, j+1] = res1

xfin[len(xfin)-1][0] = 1

mat1 = np.zeros((10000,3))
hh = 0

for j in range(len(tem)):

    xfin2 = np.zeros((len(np.transpose(xxx[:, 1])) + 1, 2))

    for i in range(len(xxx)):
        xfin2[i, 0] = xxx[i, 1]
        xfin2[i, 1] = xfin[i, j+1]

    xfin2[len(xfin2) - 1][0] = 1
    xfin2[0][1] = 0

    hull = ConvexHull(xfin2)

    for simplex in hull.simplices:

        if xfin2[simplex, 0][0] != 1 and xfin2[simplex, 0][0] != 0 and xfin2[simplex, 0][1] != 1 and xfin2[simplex, 0][
            1] != 0:

            if abs(xfin2[simplex, 0][0] - xfin2[simplex, 0][1]) > 1 / sep * 1.3:
                mat1[hh][1] = tem[j]
                mat1[hh][0] = xfin2[simplex, 0][0]
                mat1[hh][2] = xfin2[simplex, 0][1]
                hh = hh+1

f = plt.figure(figsize=(6,5))

plt.style.use('seaborn-whitegrid')

s=1.0
plt.scatter(mat1[:,0],mat1[:,1], s)

plt.xlabel('Zn Concentration (at.)')
plt.ylabel('Temperature K')

ax = plt.axes()
ax.set(xlim=[0,1], ylim=[298,1000])

# plt.show()|

f.savefig("plot.pdf")
