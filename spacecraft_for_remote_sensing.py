import math

''' I. Расчёт основных характеристик ДЗЗ '''

''' Физические константы '''
r_earth_km = 6371               # средний радиус Земли
mu_km3_s2 = 3.986e5             # гравитационный параметр Земли
wavelength_m = 1550e-9          # длина волны излучения (ИК, C-band)
speed_of_light_m_s = 299792458  # скорость света
k = 1.380609e-23                # постоянная Больцмана в Дж/К

''' Характеристики камеры и матрицы для ДЗЗ '''
fov_deg = 6.0                   # поле зрения камеры ДЗЗ
ccd_size = [12000, 1]           # размер ПЗС-матрицы в пикселях
pixel_size_m = 6.5e-6           # размер пикселя ПЗС
focal_length_m = 0.8            # фокусное расстояние камеры

''' Характеристики (лазерного) канала связи'''
bits_per_pixel = 12             # "радиометрическое" разрешение
channel = 4                     # количество спектральных каналов

''' Характеристики орбиты '''
n = 1                                                        # суточная кратность орбиты
cycles = 16                                                  # число витков
altitude_km = 42241.12 * (n/cycles)**(2/3) - r_earth_km      # высота солнечно-синхронной Земной орбиты

''' Определение параметров съёмки '''
eta_rad = math.radians(fov_deg / 2)                                     # надирный угол
epsilon_rad = math.acos(                                                # рабочий угол места    
    ((r_earth_km + altitude_km) / r_earth_km) * math.sin(eta_rad)    
)
lambda_rad = math.pi / 2 - eta_rad - epsilon_rad                        # центральный угол

fov_total_length_eta_km = 2 * eta_rad * altitude_km                     # ширина полосы обзора по надирному углу
fov_total_length_lambda_km = 2 * lambda_rad * r_earth_km                # ширина полосы обзора по центральному углу

ccd_x_size_mm = pixel_size_m * ccd_size[0] * 1000                       # поперечный размер матрицы (x)
ccd_y_size_mm = pixel_size_m * ccd_size[1] * 1000                       # продольный размер матрицы (y)

ifov_x_rad = ccd_x_size_mm / (1000 * focal_length_m)                    # (мгновенное) поперечное поле зрения
ifov_y_rad = ccd_y_size_mm / (1000 * focal_length_m)                    # (мгновенное) поперечное поле зрения

image_x_size_km = ifov_x_rad * altitude_km                              # поперечный размер кадра
image_y_size_km = ifov_y_rad * altitude_km                              # продольный размер кадра

gsd_x_m = 1000 * image_x_size_km / ccd_size[0]                          # разрешение камеры в поперечном направлении
gsd_y_m = 1000 * image_y_size_km / ccd_size[1]                          # разрешение камеры в продольном направлении

''' Определение характеристик передачи данных '''
period_s = 2 * math.pi * math.sqrt(                                     # период обращения КА
    ((r_earth_km + altitude_km) ** 3) / mu_km3_s2
)
v_shadow_km_s = 2 * math.pi * r_earth_km / period_s                     # скорость подспутниковой точки
freq_hz = v_shadow_km_s / image_y_size_km                               # частота обновления изображения                   
data_flow_megabit_s = (                                                 # выходной поток данных
    ccd_size[0] * ccd_size[1] * bits_per_pixel * channel * freq_hz
    ) / 1e6   
data_per_cycle_megabit = data_flow_megabit_s * period_s                 # объём информации, передаваемый за виток  
interval_s = 480                                                        # временной интервал в 8 мин
data_per_interval_megabit = data_flow_megabit_s * interval_s            # объём информации, передаваемый за интервал 8 мин    

''' Вывод данных в консоль '''
print()
print("---------------------------------------------")
print("Расчёт основных характеристик ДЗЗ")
print("---------------------------------------------")
print(f"высота орбиты КА: {altitude_km:.3f} км")
print()
print(f"надирный угол: {eta_rad:.6f} рад")
print(f"рабочий угол места: {epsilon_rad:.6f} рад")
print(f"центральный угол: {lambda_rad:.6f} рад")
print()
print(f"ширина полосы обзора по надирному углу: {fov_total_length_eta_km:.6f} км")
print(f"ширина полосы обзора по центральному углу: {fov_total_length_lambda_km:.6f} км")
print()
print(f"поперечный размер матрицы (x): {ccd_x_size_mm:.3e} мм")
print(f"продольный размер матрицы (y): {ccd_y_size_mm:.3e} мм")
print()
print(f"мгновенное поле зрения по x: {ifov_x_rad:.3e} рад")
print(f"мгновенное поле зрения по y: {ifov_y_rad:.3e} рад")
print()
print(f"размер кадра по x: {image_x_size_km:.6e} км")
print(f"размер кадра по y: {image_y_size_km:.6e} км")
print()
print(f"разрешение на местности (GSD) по x: {gsd_x_m:.6e} м")
print(f"разрешение на местности (GSD) по y: {gsd_y_m:.6e} м")
print()
print(f"период обращения: {period_s:.2f} с")
print(f"скорость подспутниковой точки: {v_shadow_km_s:.3f} км/с")
print(f"частота обновления изображения: {freq_hz:.3f} Гц")
print(f"выходной поток данных: {data_flow_megabit_s:.3f} Мбит/с")
print(f"объём информации, передаваемый за виток: {data_per_cycle_megabit:.3f} Мбит")
print(f"объём информации, передаваемый за 8 мнут: {data_per_interval_megabit:.3f} Мбит")
print()


''' II. Оценка бюджета лазерной линии связи КА-Земля '''

''' Характеристики линии связи КА с наземной станцией '''
epsilon_min_grad = 5.0                                                         # минимальный угол места (видимости) КА в градусах
epsilon_min_rad = math.radians(epsilon_min_grad)                               # минимальный угол места (видимости) КА в радианах

tx_freq_Hz = speed_of_light_m_s / wavelength_m        # частота передачи данных 

space_telescope_diam_m = 5e-2                                   # диаметр телескопа на КА для приёма/передачи лазерного сигнала
gs_aperture_usage_coeff = 0.6                                   # коэффициент использования апертуры наземной приеёмной станции
space_aperture_usage_coeff = 0.6                                # коэффициент использования апертуры телескопа КА (учёт центрального отверстия и других потерь)      
space_efficiency = 0.85                                         # КПД телескопа на КА (учёт потерь в оптической системе, неидеальной отражающей поверхности и других факторов)
gs_efficiency = 0.8                                             # КПД телескопа наземной приёмной станции    
space_pointing_error_rad = 5e-6                                 # ошибка наведения лазерного терминала КА
gs_pointing_error_rad = 5e-6                                    # ошибка наведения телескопа наземной станции 
gs_telescope_diam_m = 0.5                                       # диаметр телескопа наземной станции для приёма лазерного сигнала
                                       
space_laser_tx_power_mW = 1e3                            # мощность лазерного передатчика на КА в мВт
space_laser_tx_power_dBp = 10 * math.log10(              # мощность лазерного передатчика на КА в дБм (относительно 1 мВт)
    space_laser_tx_power_mW
)   

signal_energy_J_per_bit = 9                                     # энергия сигнала на бит для обеспечения требуемого Eb/N0
noise_power_density_W_per_Hz = 1                                # спектральная плотность мощности шума
required_Eb_per_N0_dB = 10 * math.log10(                        # требуемый Eb/N0 в дБ
    signal_energy_J_per_bit / noise_power_density_W_per_Hz
)  

''' Определение бюджета лазерной линии '''

''' Определение геометрических характеристик видимости КА с наземной станции '''
eta_max_rad = math.asin(                                                            # максимальный надирный угол, при котором КА всё ещё виден с наземной станции
    (r_earth_km / (r_earth_km + altitude_km)) * math.cos(epsilon_min_rad)
)
lambda_max_rad = math.pi / 2 - eta_max_rad - epsilon_min_rad                        # максимальный центральный угол, при котором КА всё ещё виден с наземной станции
distance_max_km = r_earth_km * (math.sin(lambda_max_rad) / math.sin(eta_max_rad))   # максимальное расстояние от КА до наземной станции

''' Потери в атмосфере, при лёгком дожде на трасса и в оптическом тракте '''
atmosphere_loss_dB = 0.25 * distance_max_km                     # потери в атмосфере (эмпирическая формула, зависящая от высоты)
light_rain_loss_dB = 6.5 * distance_max_km                      # потери при умеренном дожде (эмпирическая формула, зависящая от высоты)    
optical_channel_loss_dB = 5                                     # потери в оптическом тракте 

''' Расчёт характеристик телескопа на КА '''
space_telescope_geom_area_m2 = math.pi * (space_telescope_diam_m ** 2)/ 4                   # геометрическая площадь телескопа на КА                                  
space_telescope_eff_area_m2 = space_aperture_usage_coeff * space_telescope_geom_area_m2     # эффективная площадь телескопа на КА с учётом коэффициента использования апертуры
space_telescope_eff_diam_m = math.sqrt(4 * space_telescope_eff_area_m2 / math.pi)           # эффективный диаметр телескопа на КА

space_telescope_directivity = space_telescope_eff_area_m2 * 4 * math.pi * (                 # направленность телескопа на КА 
    tx_freq_Hz / speed_of_light_m_s
) ** 2                
space_telescope_directivity_dB = 10 * math.log10(space_telescope_directivity)               # направленность телескопа на КА в дБ

space_telescope_gain = space_efficiency * space_telescope_directivity                       # усиление телескопа на КА                                         
space_telescope_gain_dB = 10 * math.log10(space_telescope_gain)                             # усиление телескопа на КА в дБ

''' Расчёт характеристик телескопа на наземной станции '''
gs_telescope_geom_area_m2 = math.pi * (gs_telescope_diam_m ** 2) / 4                        # геометрическая площадь телескопа наземной станции
gs_telescope_eff_area_m2 = gs_aperture_usage_coeff * gs_telescope_geom_area_m2              # эффективная площадь телескопа наземной станции с учётом коэффициента использования апертуры
gs_telescope_eff_diam_m = math.sqrt(4 * gs_telescope_eff_area_m2 / math.pi)                 # эффективный диаметр телескопа наземной станции

gs_telescope_directivity = gs_telescope_eff_area_m2 * 4 * math.pi * (                       # направленность телескопа наземной станции
    tx_freq_Hz / speed_of_light_m_s
) ** 2        
gs_telescope_directivity_dB = 10 * math.log10(gs_telescope_directivity)                     # направленность телескопа наземной станции в дБ 

gs_telescope_gain = gs_efficiency * gs_telescope_directivity                                # усиление телескопа наземной станции
gs_telescope_gain_dB = 10 * math.log10(gs_telescope_gain)                                   # усиление телескопа наземной станции в дБ

''' Расчёт потерь на ошибку наведения лазерного терминала КА '''
theta_spacecraft_rad = wavelength_m / space_telescope_eff_diam_m                            # ширина диаграммы направленности лазерного терминала КА
space_pointing_loss_dB = 12 * (space_pointing_error_rad / theta_spacecraft_rad) ** 2        # потери на ошибку наведения лазерного терминала КА

''' Расчёт потерь на ошибку наведения телескопа наземной станции '''
theta_ground_rad = wavelength_m / gs_telescope_eff_diam_m                                   # ширина диаграммы направленности телескопа наземной станции
gs_pointing_loss_dB = 12 * (gs_pointing_error_rad / theta_ground_rad) ** 2                  # потери на ошибку наведения телескопа наземной станции 

''' Расчёт сигнала и потерь на трассе '''
loss_in_free_space_dB = 20 * math.log10(distance_max_km * 1000 * 4 * math.pi / wavelength_m)        # потери в свободном пространстве в дБ

rx_signal_power_dBm = (space_laser_tx_power_dBp + space_telescope_gain_dB + gs_telescope_gain_dB    # мощность сигнала на входе приёмника в дБм с учётом всех усилений и потерь
            - space_pointing_loss_dB - gs_pointing_loss_dB -
            - atmosphere_loss_dB - light_rain_loss_dB - optical_channel_loss_dB
            - loss_in_free_space_dB)

T_na_K = 15 + 30 / gs_telescope_eff_diam_m + 180 / epsilon_min_grad                                 # шумовая температура приёмника в К (эмпирическая формула, зависящая от эффективного диаметра приёмного телескопа и угла места)
T_na_dBK = 10 * math.log10(T_na_K)                                                                  # шумовая температура приёмника в дБК

k_dBm_per_Hz_K = 10 * math.log10(k * 1000)                                                          # постоянная Больцмана в дБм/Гц/К (перевод в мВт)   

data_flow_dB = 10 * math.log10(data_flow_megabit_s)                                                 # частота передачи данных в дБ

required_rx_signal_power_dBm = k_dBm_per_Hz_K + T_na_dBK + data_flow_dB + required_Eb_per_N0_dB     # требуемый уровень сигнала на входе приёмника в дБм для обеспечения требуемого Eb/N0 с учётом шумовой температуры и частоты передачи данных

margin_dB = rx_signal_power_dBm - required_rx_signal_power_dBm                                      # бюджет линии (запас) в дБ

print("---------------------------------------------")
print("Бюджет лазерной линии связи")
print("---------------------------------------------")
print(f": {space_laser_tx_power_dBp:.3f} дБм")
print(f"минимальный угол места: {math.degrees(epsilon_min_rad):.2f}°")
print(f"максимальная дальность: {distance_max_km:.3f} км")
print(f"частота передачи: {tx_freq_Hz:.3e} Гц")
print()
print("--- Телескоп КА ---")
print(f"диаметр телескопа КА: {space_telescope_diam_m*1000:.1f} мм")
print(f"геометрическая площадь космического телескопа: {space_telescope_geom_area_m2:.3f} м²")
print(f"эффективная площадь телескопа КА: {space_telescope_eff_area_m2:.3e} м²")
print(f"эффективный диаметр телескопа КА: {space_telescope_eff_diam_m*1000:.1f} мм")
print(f"направленность телескопа КА: {space_telescope_directivity:.3e} (раз), или {space_telescope_directivity_dB:.3f} дБ")
print(f"усиление телескопа КА: {space_telescope_gain} (раз), или {space_telescope_gain_dB:.3f} дБ")
print(f"потери на ошибку наведения КА: {space_pointing_loss_dB:.3f} дБ")
print()
print("--- Телескоп наземной станции ---")
print(f"диаметр наземного телескопа: {gs_telescope_diam_m:.2f} м")
print(f"геометрическая площадь наземного телескопа: {gs_telescope_geom_area_m2:.3f} м²")
print(f"эффективная площадь наземного телескопа: {gs_telescope_eff_area_m2:.3f} м²")
print(f"эффективный диаметр наземного телескопа: {gs_telescope_eff_diam_m:.3f} м")
print(f"направленность наземного телескопа: {gs_telescope_directivity:.3e} (раз)")
print(f"усиление наземного телескопа: {gs_telescope_gain_dB:.3f} дБ")
print(f"потери на ошибку наведения НС: {gs_pointing_loss_dB:.3f} дБ")
print()
print("--- Потери на трассе ---")
print(f"потери в свободном пространстве: {loss_in_free_space_dB:.3f} дБ")
print(f"потери в атмосфере: {atmosphere_loss_dB:.3f} дБ")
print(f"потери при лёгком дожде: {light_rain_loss_dB:.3f} дБ")
print(f"потери в оптическом тракте: {optical_channel_loss_dB:.3f} дБ")
print()
print("--- Мощность сигнала ---")
print(f"суммарная мощность на входе приёмника: {rx_signal_power_dBm:.3f} дБм")
print()
print("--- Требуемый уровень сигнала ---")
print(f"шумовая температура приёмника: {T_na_K:.3f} К ({T_na_dBK:.3f} дБК)")
print(f"постоянная Больцмана: {k_dBm_per_Hz_K:.3f} дБм/Гц/К")
print(f"требуемое Eb/N0: {required_Eb_per_N0_dB:.1f} дБ")
print(f"требуемый уровень сигнала на входе: {required_rx_signal_power_dBm:.3f} дБм")
print()
print(f"бюджет линии (запас): {margin_dB:.3f} дБ")
if margin_dB > 0:
    print("✅ запас положительный, связь возможна")
else:
    print("❌ запас отрицательный, требуется увеличение мощности или уменьшение потерь")