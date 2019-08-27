# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 19:37:11 2019

@author: aniba
"""

from matplotlib.pylab import *
import csv

a = 1.   # m, Largo del dominio
b = 0.5  # m, Alto del dominio
c = 0.6 # m, Ancho del dominio

Nx = 20 # Numero de intervalos en x, @5 cm
Ny = 10 # Numero de intervalos en y, @10 cm
Nz = 12 # Numero de intervalos en z, @5 cm

#los puntos en el espacio se encuentran a 5 cm de distancia (en multiples planos 2-D)
# IMPORTANTE = Son distintos para que los deltas sean iguales.

# Discretizacion espacial
dx = a/Nx
dy = b/Ny
dz = c/Nz

h = dx    # = dy

if dx != dy:
    print("ERRROR!!!!! dx != dy")
    exit(-1)   #-1 le dice al SO que el programa fallo.....

#Funcion de conveniencia para calcular coordenadas del punto (i,j)
coords = lambda i, j , z: (dx*i, dy*j, dz*z)
x, y, z = coords(4,2,2) 

print ("x = ", x)
print ("y = ", y)
print ("z = ", z)

u_k = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
u_km1 = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)

#CB esencial toda la matriz llena de 20
u_k[0,:,:] = 20.
u_k[:,0,:] = 20.
u_k[:,:,0] = 20.
u_k[Nx,:,:] = 20.
u_k[:,Ny,:] = 20.
u_k[:,:,Nz] = 20.

# Propiedades del hormigon
K = 116.     # m^2 / s, conductividad termica
c = 390.     # J / kg*C, calor especifico
rho = 7140.  # kg/m^3, densidad

#Arreglo de parametros
alpha_0 = 0.0001
dt = alpha_0*(c*rho*dx**2)/K
alpha = K*dt/(c*rho*dx**2)
dt = 1800 # s

#Lista de temperatura esta medida cada 80.86642599277978 segundos.

with open('TemperaturaAmbiente.csv') as csv_file:
	csv_reader = csv.reader(csv_file, delimiter = ';')
	line_count = 0
	contador = 0
	temp = []
	for row in csv_reader:
		contador += 1
		if contador < 30:
			continue
		else:
			contador = 0
			if line_count == 0:
				line_count += 1
			else:
				if len(temp) == int32(3600*24*7/dt):
					break
				temp.append(float(row[3]))
				line_count += 1


#Creando listas que guardaran los datos de las temp de cada sensor
s1 = []
s2 = []
s3 = []
s4 = []
s5 = []
s6 = []
s7 = []
s8 = []
s9 = []

temp2 = [] #temp generado en los tiempos entregados en los datos

#Loop en el tiempo 
tiempos = []
dnext_t = 18000/dt # s
next_t = 0 

# Parametros para la generacion de calor
beta = 1.05 # Pendiente de hidratacion
thau = 10.3 # Parametro de hidratacion
Cc = 160    # kg, cantidad de cemento
Tc = 20     # C, Temperatura del concreto
Tr = 23     # C, Temperatura de referencia
Cr = 23     # C,
R = 8.31    # J / mol*k, Ctte universal de los gases
E = 38380   # J / mol, Energia de activacion
H = 355.24  #J / kg, Potencial heat generation

for tiempo in range(int32(7*3600*24/dt)): #rango de 7 dias de simulacion en segundos
    t = dt*tiempo + 1
    #print ("k = ", tiempo, " t = ", t)
    Temperatura_Ambiente= temp[tiempo] #Ecuacion para la temperatura ambiente, expresado en videos de profesor
    u_k[:, :, :] = H*Cc*((thau/t)**beta)*exp(-(thau/t)**beta)*exp((E/R)*(1/(273 + Tr) - 1/(273 + Tc)))*3
    u_k[:, Ny, :] = temp[tiempo%(dt)]

    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            for r in range(1,Nz-1):
            #Laplaciano
                nabla_u_k = (u_k[i-1,j,r] + u_k[i+1,j,r] + u_k[i,j-1,r] + u_k[i,j+1,r]+ u_k[i,j,r-1]+u_k[i,j,r+1] - 6*u_k[i,j,r])/h**2
            #Forward euler..
                u_km1[i,j,r] = u_k[i,j,r] + alpha*nabla_u_k

    #CB esencial
    u_km1[0,:,:] = u_km1[1,:,:]
    u_km1[:,0,:] = u_km1[:,0,:]
    u_km1[:,:,0] = u_km1[:,:,1]
    u_km1[Nx,:,:] = u_km1[Nx-1,:,:]
    u_km1[:,:,Nz] = u_km1[:,:,Nz-1]
    #Avanzar la solucion a k + 1
    u_k = u_km1

    #Temperatura ambiente se agrega en parte superior de bloque, en punto maximo de y

    if tiempo >= next_t:
        next_t += dnext_t
        s1.append(u_k[10,0,0])
        s2.append(u_k[10,5,0])
        s3.append(u_k[10,10,0])
        s4.append(u_k[0,0,0])
        s5.append(u_k[0,5,0])
        s6.append(u_k[0,10,0])
        s7.append(u_k[10,0,6])
        s8.append(u_k[10,5,6])
        s9.append(u_k[10,10,6])
        temp2.append(temp[tiempo%dt])
        tiempos.append(t)


#Ploteo de los puntos 
#Ploteo de los puntos 
xlabel('Tiempo (s)')
ylabel('Temperatura (C)')
title('Semana 1')
plot(tiempos, s1, label = 's1')
plot(tiempos, s2, label = 's2')
plot(tiempos, s3, label = 's3')
plot(tiempos, s4, label = 's4')
plot(tiempos, s5, label = 's5')
plot(tiempos, s6, label = 's6')
plot(tiempos, s7, label = 's7')
plot(tiempos, s8, label = 's8')
plot(tiempos, s9, label = 's9')
plot(tiempos, temp2, label = 'Ambiente')
legend()
savefig('grafico.png')