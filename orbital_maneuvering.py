import math

r_earth_km = 6371               # средний радиус Земли
mu_km3_s2 = 3.986e5             # гравитационный параметр Земли

target_altitude_km = 282          # целевая орбита (высота над поверхностью Земли)
launch_orbit_km = 400             # начальная орбита (высота над поверхностью Земли)

r_1_km = r_earth_km + target_altitude_km         # радиус целевой орбиты
r_2_km = r_earth_km + launch_orbit_km            # радиус начальной орбиты 

vel_orb_1_km_s = math.sqrt(mu_km3_s2 / r_1_km)       # орбитальная скорость на целевой орбите
vel_orb_2_km_s = math.sqrt(mu_km3_s2 / r_2_km)       # орбитальная скорость на начальной орбите

vel_perigee_km_s = vel_orb_1_km_s * math.sqrt(2 * r_2_km / (r_1_km + r_2_km))  # скорость в перигее эллиптической орбиты
vel_apogee_km_s = vel_orb_2_km_s * math.sqrt(2 * r_1_km / (r_1_km + r_2_km))    # скорость в апогее эллиптической орбиты

delta_v1_km_s = vel_perigee_km_s - vel_orb_1_km_s  # ΔV для перехода на низкую круговую орбиту
delta_v2_km_s = vel_orb_2_km_s - vel_apogee_km_s  # ΔV для перехода на эллиптическую орбиту

delta_v_total = delta_v1_km_s + delta_v2_km_s

p = math.pi * (0.5 * (r_1_km + r_2_km)) ** (3 / 2) / math.sqrt(mu_km3_s2)  # время перехода

print("=== Орбитальный маневр: Hohmann Transfer ===")
print(f"Целевая орбита: {target_altitude_km} км над поверхностью Земли")
print(f"Начальная орбита: {launch_orbit_km} км над поверхностью Земли")
print(f"ΔV для перехода на низкую круговую орбиту: {delta_v1_km_s:.3f} км/с")       
print(f"ΔV для перехода на эллиптическую орбиту: {delta_v2_km_s:.3f} км/с")
print(f"Общий ΔV: {delta_v_total:.3f} км/с")
print(f"Время перехода (полпериода эллиптической орбиты): {p:.2f} секунд ({p/60:.2f} минут)")