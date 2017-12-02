import numpy as np
import matplotlib.pyplot as plt
from pylab import *

datos=np.genfromtxt('RadialVelocities.dat').T
datos2=np.genfromtxt('daticos.txt').T



plt.plot(datos[0],datos[1],color='blue',label='Velocidades observadas')
plt.plot(datos[0],datos2,color='red',label='Prediccion')
plt.title('Velocidad vs radio')
plt.xlabel('Radio [kpc]')
plt.ylabel('Velocidad [km/s]')
plt.legend()
plt.savefig('curvaRotacion.png')

