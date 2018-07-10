import numpy as np
import matplotlib.pyplot as plt
data=np.loadtxt("1.txt")
data1=np.loadtxt("2.txt")
plt.figure()
x=[]
y=[]
fig = plt.figure()
for i in range(15):
        for j in range(15):
                if (data[i][j]==1.0):
                        x.append(i)
                        y.append(j)

plt.plot(data1[:,0],data1[:,1],'.-')
plt.hold('on')
plt.grid('on')
plt.scatter(x,y)
plt.hold('on')
plt.grid('on')
plt.show()