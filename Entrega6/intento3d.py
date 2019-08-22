from matplotlib.pylab import *


a = 10.    		#Ancho del dominio
b = 10.
c = 10.
          #Largo del dominio
Nx = 40   		#Numero de intervalos en x
Ny = 40    		#Numero de intervalos en Y
Nz = 40
#Valores para los cuales no se malogra el grafico de la disipacion de temperatura

dx = b / Nx   	#Discretizacion espacial en X
dy = a / Ny   	#Discretizacion espacial en Y
dz = c/ Nz

h = dx    # = dy


if dx != dy:
	print("ERRROR!!!!! dx != dy")
	exit(-1)   #-1 le dice al SO que el programa fallo.....

#Funcion de conveniencia para calcular coordenadas del punto (i,j)


coords = lambda i, j , z: (dx*i, dy*j, dz*z)
x, y , z = coords(4,2,2) 

print "x = ", x
print "y = ", y
print "z = ", z

u_k = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
u_km1 = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)

#CB esencial
u_k[0,:,:] = 20.
#u_k[:,0,:] = 20.
u_k[:,:,0] = 20.
u_k[Nx,:,:] = 20.
u_k[:,Ny,:] = 20.
u_k[:,:,Nz] = 20.


#Primeras 3 condiciones de borde son replicadas cosntantemente a lo largo del codigo
#Ultima condicion es extraida luego de que empieza a aplicarse temperatura ambiente

##Buena idea definir funciones que hagan el codigo expresivo
def printbien(u):
	print u.T[Nx::-1,:]
            #Imprime con el eje y invertido
printbien(u_k)

def imshowbien(u):
	imshow(u.T[Nx::-1,:,:],vmin=0,vmax=40) #vmin y vmax definen los valores extremos para la temperatura

#Parametros del problema (hierro)
dt = 60.       # s # tiempo 1 minuto
K = 79.5       # m^2 / s   
c = 450.       # J / kg C
rho = 7800.    # kg / m^3
alpha = K*dt/(c*rho*dx**2)

# dx =  0.166666666667
# dt = 1.0
# alpha =  0.000815384615385

alpha_bueno = 0.0001
#dt = alpha_bueno*(c*rho*dx**2)/K # no se utiliza, dt=1 minuto
alpha = K*dt/(c*rho*dx**2)
Tiempo=24*60*60 #1 dia en segundos
#Temperatura_Ambiente=20+10*sin((2*pi/Tiempo)*t) #Ecuacion para la temperatura ambiente

#Informar cosas interesantes
print "dt = ", dt
print "dx = ", dx
print "K = ", K
print "c = ", c
print "rho = ", rho
print "alpha = ", alpha

k = 0

# figure(1)
# imshowbien(u_k)
# title("k = {}   t = {} s".format(k, k*dt))
# savefig("movie/frame_{0:04.0f}.png".format(k))
# close(1)

#Loop en el tiempo 

for k in range(5000): #rango arbitrario de prueba
	t = dt*(k+1)
	print "k = ", k, " t = ", t
	Temperatura_Ambiente=20+10*sin((2*pi/Tiempo)*t) #Ecuacion para la temperatura ambiente, expresado en videos de profesor
	#CB esencial
	#CB esencial
	u_k[0,:,:] = 20.
	#u_k[:,0,:] = 20.
	u_k[:,:,0] = 20.
	u_k[Nx,:,:] = 20.
	u_k[:,Ny,:] = 20.
	u_k[:,:,Nz] = 20.
	#Ultima condicion de borde no es replicada despues de agregar calor
	u_k[:,0,:] = Temperatura_Ambiente
	#Temperatura ambiente se agrega en parte superior de bloque, en punto maximo de y
	#Loop en el espacio   i = 1 ... n-1   u_km1[0] = 0  u_km1[n] = 20
	for i in range(1,Nx):
		for j in range(1,Ny):
			for r in range(1,Nz):
			#Algoritmo de diferencias finitas 2-D para difusion

			#Laplaciano
				nabla_u_k = (u_k[i-1,j,r] + u_k[i+1,j,r] + u_k[i,j-1,r] + u_k[i,j+1,r]+ u_k[i,j,r-1]+u_k[i,j,r+1] -6*u_k[i,j,r])/h**2

			#Forward euler..
				u_km1[i,j,r] = u_k[i,j,r] + alpha*nabla_u_k

	#hay que revisar lo que sigue...

	#CB esencial
	u_km1[0,:,:] = u_km1[1,:,:]
	u_km1[Nx,:,:] = u_km1[Nx-1,:,:]
	#u_k[:,0,:] = 20.
	u_km1[:,:,0] = u_km1[:,:,1]
	u_km1[:,:,Nz] = u_km1[:,:,Nz-1]
	u_km1[:,Ny,:] = u_km1[:,Ny-1,:]
	#Avanzar la solucion a k + 1
	u_k = u_km1

	#CB esencial una ultima vez
	u_k[0,:] = 20.
	u_k[:,0] = 20.
	u_k[Nx,:] = 20.
	#Ultima condicion de borde no es replicada despues de agregar calor
	u_k[:,0,:] = Temperatura_Ambiente
	#Temperatura ambiente se agrega en parte superior de bloque, en punto maximo de y

	#print "Tmax = ", u_k.max()7



	#figure(2)
	#imshowbien(u_k)
	#title("k = {}   t = {} s".format(k, (k+1)*dt))
	#colorbar()
	#if k%20==0: # Se guarda imagen cada 20 frames
	#	savefig("movie1/frame_{0:05.0f}.png".format(k))
	close()

