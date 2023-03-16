import matplotlib.pyplot as plt
import numpy as np

p = [1,2,4,8,16]
t = [[0.039  , 0.040 , 0.044 , 0.039 ,0.042],
     [0.045  , 0.041 , 0.034 , 0.030 , 0.027],
     [0.079  , 0.061 , 0.048 , 0.041 , 0.042],
     [46.57  , 36.48 , 26.28 , 17.48 , 11.13],
     [0.0017  , 0.0022 , 0.0030 , 0.0038 , 0.0042],
     [1.241  , 0.747 , 0.491 , 0.333 , 0.294],
     [0.91  , 1.04 , 1.03 , 0.85 , 0.71],
     [2.80  , 3.44 , 3.81 , 3.46 , 2.99],
     [3.40  , 3.00 , 2.47 , 1.84 , 1.49]]
s = []
e = []

for i in range(9):
    temp = [0,0,0,0,0]
    for j in range(5):
        temp[j] = t[i][0] / t[i][j]
    s.append(temp)
for i in range(9):
    temp = [0,0,0,0,0]
    for j in range(5):
        temp[j] = s[i][j] / p[j]
    e.append(temp)

def printl(l):
    for i in range(len(l)):
        if(i == len(l)-1):
            print("%.3f"%l[i],end=" \\\\")
        else:
            print("%.3f"%l[i],end=" & ")

for i in range(9):
    print("\\hline")
    print("Test",end="")
    print(i,end=" & ")
    printl(e[i])
    print("\n")
print("\\hline\n")
# xpoints = np.array(p)
# ypoints = np.array(e1)
# zpoints = np.array(e2)
# wpoints = np.array(e3)
# vpoints = np.array(e5)
# upoints = np.array(e8)

# plt.plot(xpoints,ypoints, marker = '.',label= "testcase 1")
# plt.plot(xpoints,zpoints, marker = '.',label= "testcase 2")
# plt.plot(xpoints,wpoints, marker = '.',label= "testcase 3")
# plt.plot(xpoints,vpoints, marker = '.',label= "testcase 5")
# plt.plot(xpoints,upoints, marker = '.',label= "testcase 8")
# plt.xlabel("Number of Nodes")
# plt.ylabel("Efficieny")
# plt.title("Effiiency vs Nodes")
# plt.legend()
# plt.show()