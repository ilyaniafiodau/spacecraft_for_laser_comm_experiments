import numpy as np
from matplotlib import pyplot as plt

r_earth_km = 6371.0                                                        # радиус Земли в километрах
n = 1                                                                      # суточная кратность орбиты
cycles = np.arange(10, 18, 1)                                               # число витков
altitude_km = np.array(
    [42241.12 * (n/cycles)**(2/3) - r_earth_km for cycles in cycles]
)      # высота солнечно-синхронной Земной орбиты

plt.figure()
plt.xlabel("Число витков")
plt.ylabel("Высота орбиты, км")
plt.minorticks_on()
plt.xlim(10, 17)
plt.yticks(np.arange(0, 3000, 250))
plt.grid(which = 'major')
plt.grid(which = 'minor', linestyle = ':')
plt.plot(cycles, altitude_km)
plt.title("Зависимость высоты орбиты КА от числа витков из условия единичной суточной кратности ССО")
plt.show()