#Proyecto 1 MCOC
#
from numpy import *
from matplotlib.pylab import *

L = 1.		#largo del dominio
n = 100		#numero de intervalos

dx = L / n 		#discretizacion espacial

#Vector con todos los x... puntos del espacio
x = linspace(0,L,n+1)

#condicion inicial
def fun_u0(x):
	return 10*exp(-(x-0.5)**2/0.1**2) 

u0 = fun_u0(x)

k=0

#Creando el vector solucion u en el tiempo o paso k
u_k = u0.copy() #crea una nueva instancia del vector en la memoria

#Condiciones de borde (esenciales)
u_k[0]=0.
u_k[n]=20.

#Temperatura en el tiempo k + 1 = dt * (k+1) 
u_km1 = u_k.copy()

#Parametros del problema
dt = 1		#s
K = 79.5	#m2/s
rho = 7800	#J / kg*C
c = 450		#kg / m3

print "dt = ", dt
print "dx = ", dx
print "K = ", K
print "c = ", c
print "rho = ", rho

plot(x,u0, "k--")

#Loop en elt iempo
k = 0
for k in range(3600*5):
	t = dt*k
	print "k =", k, "t = ", t
	#Loop en el espacio i = 1 ... n-1   u_km1[0] = 0  u_km1[n] = 20
	u_k[0] = 0
	u_k[n] = 20
	u_km1[1] = 0 + K*dt/(c*rho*dx**2)*(u_k[1+1] - 2*u_k[1] + u_k[1-1])
	for i in range(2,n):
		#Algoritmo de diferencia finitas 1-D para difusion
		u_km1[i] = u_k[i] + K*dt/(c*rho*dx**2)*(u_k[i+1] - 2*u_k[i] + u_k[i-1])
	#Avanzar la solucion a k+1
	u_k = u_km1
	if k % 200 == 0:
		plot (x,u_k)

title("k = {}   t = {} s".format(k, k*dt))
show ()
