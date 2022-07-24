"""A code to draw geodesics of a sphere using 4-rank-runge-kutta"""

from sympy import *
import numpy as np
import matplotlib.pyplot as plt 

coor=symbols('x y')

""" X is a function from R^2 to R^3 the inverse of chart"""
X=[coor[0],coor[1],2*coor[0]*coor[1]+2*coor[1]-coor[0]**2-coor[1]**2]

""" Definition of Riemannian metric"""
def Metric(X):
    Xx=[diff(X[i],coor[0]) for i in range(3)]
    Xy=[diff(X[i],coor[1]) for i in range(3)]
    """Components of metric in the same order that they appear"""
    E=sum([Xx[i]*Xx[i] for i in range(3)])
    F=sum([Xx[i]*Xy[i] for i in range(3)])
    G=sum([Xy[i]*Xy[i] for i in range(3)])
    """Riemannianmetric"""
    g=np.array([[E,F],
               [F ,G]])
    ginv=1/(E*G-F**2)*np.array([[G, -F],[-F, E]])
    return g,ginv



"""Defining Christoffel symbols  The first component is for up index, the intermediate 
index is for rows and the last is for column"""


Christoffels=[]
for i in range(2):
    Christoffels.append(zeros(2,2))
    

for i in range(2):
    for k in range(2):
        for l in range(2):
            Christoffels[i][k,l]=1/2*sum([Metric(X)[1][i][m]*(diff(Metric(X)[0][m][k],coor[l])+diff(Metric(X)[0][m][l],coor[k])-diff(Metric(X)[0][k][l],coor[m])) for m in range(2)])





"""End of Christoffels"""


"""gama3 is derivative of gamma1, and gamma4 is derivative of gamma2 and our geodesic
is gama=(gamma1  ,  gama2)"""

"""
Initialization.
Note that gama contains the coordinates of geodesic and therefore gama=[[0],[0]] is the start point of it

"""
gama=[[0],[0]]
gama3=[0.02]
gama4=[1]


h=0.01

for i in range(100):
    gama3.append(gama3[i] + h*sum([-Christoffels[0][m,n].subs([(coor[0],gama[0][i]),(coor[1],gama[1][i])]) *gama[m][i] *gama[n][i] for m in range(2) for n in range(2)]) )
    gama4.append(gama4[i] + h*sum([-Christoffels[1][m,n].subs([(coor[0],gama[0][i]),(coor[1],gama[1][i])]) *gama[m][i] *gama[n][i] for m in range(2) for n in range(2)]))
    gama[0].append(gama[0][i]+h*gama3[i])
    gama[1].append(gama[1][i]+h*gama4[i])




"""Now I want to provide the list of points of geodesic"""


zc=[]
for i in range(len(gama[0])):
    zc.append(X[2].subs([(coor[0],gama[0][i]),(coor[1],gama[1][i])]))
    


teta=np.linspace(0,2*np.pi,1000)
phi=np.linspace(0,np.pi,1000)
x=np.outer(np.cos(teta),np.sin(phi))
y=np.outer(np.sin(teta),np.sin(phi))
z=np.outer(np.ones(np.size(teta)),np.cos(phi))


fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.plot(gama[0],gama[1],zc,color='r')

ax.plot_surface(x,y,z,color='green',alpha=0.4)
plt.show()


