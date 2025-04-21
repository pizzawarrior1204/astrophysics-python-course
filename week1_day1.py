import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(0, 2*np.pi, 100)
y=np.sin(x)


plt.plot(x, y,label='Sine Wave', color='red')
plt.title('Simpple Sine Wave')
plt.xlabel('X')
plt.ylabel('Sin(X)')
plt.legend()
plt.grid(True)
plt.show()