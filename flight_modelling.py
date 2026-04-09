import numpy as np
from matplotlib import pyplot as plt
import math

mu_m3s2 = 398600.4415e+09   
r_earth_m = 6371e+03              
h0_m = 282e+03                    
i0_deg = 95.68                     
F107 = 150          
F81 = F107
K0 = 1
K1 = 0
K2 = 0
K3 = 0
K4 = 0
c_ball_m2kg = np.arange(0.001, 0.006, 0.001)

i0_rad = math.radians(i0_deg)
r0_m = r_earth_m + h0_m        
V0_circular_ms = math.sqrt(mu_m3s2/r0_m)    
print(f'V0 circular = {V0_circular_ms}')
Vy0_ms = V0_circular_ms * math.cos(i0_rad)
Vz0_ms = V0_circular_ms * math.sin(i0_rad)

t0_s = 0.0
state_vector_0 = np.array(
    [
        r0_m, 0.0, 0.0,           
        0.0, Vy0_ms, Vz0_ms           
    ], dtype=np.float64
)

def density(h_m: np.float64) -> np.float64:
    h_km = h_m / 1000
    if h_km >= 120:
        dens0_kgm3 = 1.58868e-08
        a0 = 29.6418
        a1_km_1 = -0.514957
        a2_km_2 = 0.00341926
        a3_km_3 = -1.25785e-05
        a4_km_4 = 2.5727e-08
        a5_km_5 = -2.75874e-11
        a6_km_6 = 1.21091e-14
        dens_n_kgm3 = dens0_kgm3 * math.exp(a0 + a1_km_1*h_km + a2_km_2*h_km**2 + a3_km_3*h_km**3
                                  + a4_km_4*h_km**4 + a5_km_5*h_km**5 + a6_km_6*h_km**6)
        dens_kgm3 = dens_n_kgm3 * K0 * (1 + K1 + K2 + K3 + K4)
    else:
        a0i_kgm3 = 3.66e-07
        k1i_km_1 = -0.18553
        k2i_km_2 = 1.5397e-03
        h_i_km = 110
        dens_kgm3 = a0i_kgm3 * math.exp(k1i_km_1*(h_km-h_i_km) + k2i_km_2*(h_km-h_i_km)**2)
    return dens_kgm3

def right_sides (
        t: np.float64,                  
        state_vector: np.ndarray,
        i: np.int16      
        ) -> np.ndarray: 
    global mu_m3s2          
    r_m = np.linalg.norm(state_vector[0:3])   
    v_ms = np.linalg.norm(state_vector[3:6])   
    h_m = r_m - r_earth_m                      
    derivs = np.zeros(shape=[6])   
    derivs[0:3] = state_vector[3:6]     
    derivs[3:6] = - ((mu_m3s2 / r_m**3) * state_vector[0:3]) - (c_ball_m2kg[i] * density(h_m) * v_ms**2  * (state_vector[3:6] / v_ms))   
    return derivs

def right_sides_airless (
        t: np.float64,                  
        state_vector: np.ndarray        
        ) -> np.ndarray: 
    global mu_m3s2          
    r_m = np.linalg.norm(state_vector[0:3])   
    derivs = np.zeros(shape=[6])  
    derivs[0:3] = state_vector[3:6]    
    derivs[3:6] = - ((mu_m3s2 / r_m**3) * state_vector[0:3])  
    return derivs

# print("right sides = ", right_sides(t0_s, state_vector_0))

def RK4_airless(t: np.float64, dh: np.float64, state_vector: np.ndarray) -> np.ndarray:
    k1 = right_sides_airless(t, state_vector)                               
    k2 = right_sides_airless(t + dh/2, state_vector + dh * k1/2)        
    k3 = right_sides_airless(t + dh/2, state_vector + dh * k2/2)        
    k4 = right_sides_airless(t + dh, state_vector + dh * k3)              
    state_vector_next = state_vector + dh * (k1 + 2*k2 + 2*k3 + k4) / 6
    return state_vector_next

def RK4(t: np.float64, dh: np.float64, state_vector: np.ndarray, i: np.int16) -> np.ndarray:
    k1 = right_sides(t, state_vector, i)                               
    k2 = right_sides(t + dh/2, state_vector + dh * k1/2, i)        
    k3 = right_sides(t + dh/2, state_vector + dh * k2/2, i)        
    k4 = right_sides(t + dh, state_vector + dh * k3, i)              
    state_vector_next = state_vector + dh * (k1 + 2*k2 + 2*k3 + k4) / 6
    return state_vector_next
                    
dt_s = 1.0                                          
altitude0_km = (np.linalg.norm(state_vector_0[0:3]) - r_earth_m) / 1000 

def check_step_size(dh: np.float64):
    traj_points_check = [state_vector_0.copy(), ]
    state_vector_local = state_vector_0.copy()
    t_local_s = t0_s
    revs = 0
    while revs < 100:  
        prev_z = state_vector_local[2]  
        next_state = RK4_airless(t_local_s, dh, state_vector_local)
        next_z = next_state[2]
        if prev_z <= 0 and next_z >= 0:
            revs += 1
            print(f'Число витков: {revs}')
        state_vector_local = next_state
        t_local_s = t_local_s + dh
        traj_points_check.append(state_vector_local.copy())
    print(f'Последняя точка = {traj_points_check[-1]}')
    interp_point = traj_points_check[-2][0:3] + ((traj_points_check[-1][0:3] - traj_points_check[-2][0:3]) / (
        traj_points_check[-1][2] - traj_points_check[-2][2]) * (0 - traj_points_check[-2][2]))
    print(f'Интерполированная точка = {interp_point}')
    print(f'Расстояние между конечной и начальной точками восходящего узла, выраженное в метрах = ',
        np.linalg.norm(interp_point - state_vector_0[0:3]))
    
def math_model_decline_by_10km(i: np.int16) -> np.ndarray:
    print(f"Расчёт движения КА при снижении высоты орбиты на 10 км от начальной для баллистического коэффициента {c_ball_m2kg[i]} м^2/кг:")
    t_local_s = t0_s
    state_vector_local = state_vector_0.copy()
    traj_points_local = [state_vector_local.copy()]
    while (r0_m - np.linalg.norm(state_vector_local[0:3])) < 10e+03:  
        altitude_km = (np.linalg.norm(state_vector_local[0:3]) - r_earth_m) / 1000    
        if t_local_s % 5000 == 0:
            print (
                "t:", t_local_s, "с,",
                "x:", round(state_vector_local[0]/1000, 3), "км,",
                "y:", round(state_vector_local[1]/1000, 3), "км,",
                "z:", round(state_vector_local[2]/1000, 3), "км,",
                "Vx:", round(state_vector_local[3], 1), "м/с,",
                "Vy:", round(state_vector_local[4], 1), "м/с,",
                "Vz:", round(state_vector_local[5], 1), "м/с,",
                "h:", round(altitude_km, 3), "км;"
            )   
        state_vector_local = RK4(t_local_s, dt_s, state_vector_local, i)   
        t_local_s = t_local_s + dt_s
        traj_points_local.append(state_vector_local.copy())
    print (
        "t:", t_local_s, "с,",
        "x:", round(state_vector_local[0]/1000, 3), "км,",
        "y:", round(state_vector_local[1]/1000, 3), "км,",
        "z:", round(state_vector_local[2]/1000, 3), "км,",
        "Vx:", round(state_vector_local[3], 1), "м/с,",
        "Vy:", round(state_vector_local[4], 1), "м/с,",
        "Vz:", round(state_vector_local[5], 1), "м/с,",
        "h:", round(altitude_km, 3), "км;"
    )
    print("Расчёт завершён: снижение высоты орбиты на 10 км от начальной")
    return np.array(traj_points_local, dtype=np.float64)

def math_model_flight_over_100km(i: np.int16):
    print(f"Расчёт движения КА на высотах свыше 100 км для баллистического коэффициента {c_ball_m2kg[i]}:")
    t_local_s = t0_s
    state_vector_local = state_vector_0.copy()
    times_local = [t_local_s]
    heights_local = [(np.linalg.norm(state_vector_local[0:3]) - r_earth_m) / 1000]
    while (np.linalg.norm(state_vector_local[0:3]) - r_earth_m) >= 100e+03:  
        altitude_km = (np.linalg.norm(state_vector_local[0:3]) - r_earth_m) / 1000 
        if t_local_s % 15000 == 0:
            print ("t:", t_local_s, "с,",
                "x:", round(state_vector_local[0]/1000, 3), "км,",
                "y:", round(state_vector_local[1]/1000, 3), "км,",
                "z:", round(state_vector_local[2]/1000, 3), "км,",
                "Vx:", round(state_vector_local[3], 1), "м/с,",
                "Vy:", round(state_vector_local[4], 1), "м/с,",
                "Vz:", round(state_vector_local[5], 1), "м/с,",
                "h:", round(altitude_km, 3), "км;")
        state_vector_local = RK4(t_local_s, dt_s, state_vector_local, i)
        t_local_s = t_local_s + dt_s
        heights_local.append(altitude_km)
        times_local.append(t_local_s)
    print ("t:", t_local_s, "с,",
                "x:", round(state_vector_local[0]/1000, 3), "км,",
                "y:", round(state_vector_local[1]/1000, 3), "км,",
                "z:", round(state_vector_local[2]/1000, 3), "км,",
                "Vx:", round(state_vector_local[3], 1), "м/с,",
                "Vy:", round(state_vector_local[4], 1), "м/с,",
                "Vz:", round(state_vector_local[5], 1), "м/с,",
                "h:", round(altitude_km, 3), "км;")
    print (f'Время существования КА на высотах выше 100 км: {t_local_s}') 
    return np.array(times_local, dtype=np.float64), np.array(heights_local, dtype=np.float64)

check_step_size(dt_s)
traj_points_np = []
times_np = []
heights_np = []

for i in range(c_ball_m2kg.size):
   # traj_points_np.append(math_model_decline_by_10km(i))
    times_i, heights_i = math_model_flight_over_100km(i)
    times_np.append(times_i / 3600 / 24)
    heights_np.append(heights_i)

plt.figure()
plt.xlabel("Время, сутки")
plt.ylabel("Высота орбиты, км")
plt.minorticks_on()
plt.xlim([0., 140.])
plt.ylim([100., 300.])
plt.grid(which = 'major')
plt.grid(which = 'minor', linestyle = ':')
for i in range(c_ball_m2kg.size):
    plt.plot(times_np[i], heights_np[i], label = f"Баллистический коэффициент = {c_ball_m2kg[i]} м²/кг")
plt.axhline(y=272.0, color='k', linestyle='--', linewidth=1.0, label="Высота в 272 км")
plt.legend()
plt.title("Зависимость высоты орбиты от времени для различных баллистических коэффициентов")

plt.figure()
plt.xlabel("Время, сутки")
plt.ylabel("Высота орбиты, км")
plt.minorticks_on()
plt.xlim([0., 70.])
plt.ylim([260., 290.])
plt.grid(which = 'major')
plt.grid(which = 'minor', linestyle = ':')
for i in range(c_ball_m2kg.size):
    plt.plot(times_np[i], heights_np[i], label = f"Баллистический коэффициент = {c_ball_m2kg[i]} м²/кг")
plt.axhline(y=272.0, color='k', linestyle='--', linewidth=1.0, label="Высота в 272 км")
plt.legend()
plt.title("Зависимость высоты орбиты от времени для различных баллистических коэффициентов")

# ax = plt.figure().add_subplot(111, projection='3d')
# plt.xlabel("x, km")
# plt.ylabel("y, km")
# ax.set_zlabel("z, km")
# for i in range(c_ball_m2kg.size):
#     all_coordinates_km = traj_points_np[i][:, 0:3] / 1000
#     ax.plot(
#         all_coordinates_km[:, 0],
#         all_coordinates_km[:, 1], 
#         all_coordinates_km[:, 2],
#         label = f"Баллистический коэффициент = {c_ball_m2kg[i]} м^2/кг"
#     )   
# ax.legend()   
# plt.axis('equal')
# plt.grid()
# plt.title("Траектории движения КА для различных баллистических коэффициентов")

plt.show()