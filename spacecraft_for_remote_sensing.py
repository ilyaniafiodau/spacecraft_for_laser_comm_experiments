import math

''' Физические константы и исходные параметры '''
r_earth_km = 6371               # средний радиус Земли
mu_km3_s2 = 3.986e5             # гравитационный параметр Земли
wavelength_m = 1550e-9          # длина волны излучения (ИК, C-band)
speed_of_light_m_s = 299792458  # скорость света

''' I. Расчёт основных характеристик ДЗЗ '''

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
interval_s = 373.248                                                    # временной интервал 
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
print(f"объём информации, передаваемый за {interval_s} c: {data_per_interval_megabit:.3f} Мбит")
print()


# =============================================================================
# II. Оценка бюджета лазерной линии связи КА–Земля
#
# Методология: Giggenbach et al., "Link Budget Calculation in Optical LEO
# Satellite Downlinks with On/Off-Keying and Large Signal Divergence"
#
# Суммарный бюджет (формула 1 из статьи):
#   p_Rx = p_Tx + a_Tx + g_Tx + a_BW + a_FSL + a_Atm + a_Sci + g_Rx + a_Rx
#
# Усиление передающего телескопа определяется через угол расходимости пучка
# (формула 10), а не через площадь апертуры и частоту (это формула для СВЧ).
# Чувствительность приёмника задаётся числом фотонов на бит (формула 20),
# а не через шумовую температуру (это подход для СВЧ).
# =============================================================================
 
print("---------------------------------------------")
print("Бюджет лазерной линии связи КА–Земля")
print("---------------------------------------------")
 
# ── Исходные параметры линии связи ───────────────────────────────────────────
 
epsilon_min_grad = 5.0                              # минимальный угол места (°)
epsilon_min_rad  = math.radians(epsilon_min_grad)
 
space_laser_tx_power_W   = 1.0                      # мощность лазерного передатчика на КА, Вт
space_laser_tx_power_mW  = space_laser_tx_power_W * 1e3
space_laser_tx_power_dBm = 10 * math.log10(space_laser_tx_power_mW)
 
# Диаметры телескопов
space_telescope_diam_m = 5e-2                       # диаметр передающего телескопа КА, м (50 мм)
gs_telescope_diam_m    = 0.5                        # диаметр приёмного телескопа наземной станции, м (500 мм)
 
# КПД оптических трактов
space_efficiency = 0.85                             # КПД передающего оптического тракта КА (aTx в линейном виде)
gs_efficiency    = 0.80                             # КПД приёмного оптического тракта НС    (aRx в линейном виде)
 
# Коэффициент использования апертуры (учёт центрального экранирования и обрезки пучка)
space_aperture_usage_coeff = 0.60
gs_aperture_usage_coeff    = 0.60
 
# Ошибки наведения (среднеквадратические, в радианах)
space_pointing_error_sigma_rad = 5e-6               # σ_BW для лазерного терминала КА
gs_pointing_error_sigma_rad    = 5e-6               # σ для наземной станции (учитывается отдельно)
 
# Параметры чувствительности приёмника (Giggenbach, формула 20 и Таблица IV)
photons_per_bit = 250                               # типовое значение для оптимизированного InGaAs-APD при BER=1e-3
 
# Параметр атмосферного ослабления: зенитный коэффициент пропускания Tz
# Из Таблицы II статьи Giggenbach: λ=1550 нм, уровень моря, тропики/городская застройка → Tz=0.891
# Это наихудший сценарий (сценарий B). Для высокогорной чистой атмосферы (сценарий A): Tz=0.986
Tz = 0.891                                          # зенитный коэффициент пропускания атмосферы (линейный)
 
# Потери на турбулентные замирания (aSci) — принимаются равными нулю при усреднении
# за 100 мс (единичное среднее), как в примере Giggenbach (Таблица V)
a_sci_dB = 0.0
 
# Потери в оптическом тракте наведения и расщепителях внутри НС (aRx, кроме эффективности телескопа)
a_rx_internal_dB = -4.1                             # по аналогии с Таблицей V (SOFA)
 
# ── 1. Геометрия: максимальная дальность при минимальном угле места ──────────
#
# Из треугольника КА–НС–центр Земли (Giggenbach, формула 14, упрощение при HGS≈0):
#   L = sqrt((RE·sin(ε))² + 2·H0·RE + H0²) − RE·sin(ε)
# где H0 = altitude_km, RE = r_earth_km, ε = epsilon_min_rad
 
H0 = altitude_km
RE = r_earth_km
 
L_max_km = (
    math.sqrt((RE * math.sin(epsilon_min_rad))**2 + 2 * H0 * RE + H0**2)
    - RE * math.sin(epsilon_min_rad)
)
L_max_m = L_max_km * 1e3
 
print(f"минимальный угол места: {epsilon_min_grad:.1f}°")
print(f"максимальная дальность (при ε_min): {L_max_km:.3f} км")
print()
 
# ── 2. Потери в свободном пространстве — формула (13) Giggenbach ─────────────
#
#   a_FSL = 10·log10( (λ / (4π·L))² )
#         = 20·log10( λ / (4π·L) )       [отрицательное значение]
 
a_FSL_dB = 20 * math.log10(wavelength_m / (4 * math.pi * L_max_m))
 
print(f"--- Потери в свободном пространстве ---")
print(f"a_FSL = {a_FSL_dB:.3f} дБ")
print()
 
# ── 3. Атмосферное ослабление — формула (16) Giggenbach ─────────────────────
#
#   a_Atm = 10·log10( Tz^(1/sin(ε)) )
# При ε = epsilon_min_rad это наихудший (наибольший) случай ослабления.
 
a_Atm_dB = 10 * math.log10(Tz ** (1.0 / math.sin(epsilon_min_rad)))
 
print(f"--- Атмосферное ослабление ---")
print(f"зенитный коэффициент пропускания Tz: {Tz}")
print(f"a_Atm при ε={epsilon_min_grad}°: {a_Atm_dB:.3f} дБ")
print()
 
# ── 4. Передающий телескоп КА ────────────────────────────────────────────────
#
# Угол расходимости гауссова пучка по уровню 1/e² (формула 8 Giggenbach):
#   θ_e-2 = 2λ / (π · ω0)
# где ω0 = D_e-2 / 2 — радиус пучка на выходе.
#
# Соотношение апертуры и диаметра пучка (раздел III.A Giggenbach):
#   D_Tx = √2 · D_e-2  →  D_e-2 = D_Tx / √2
#   ω0 = D_e-2 / 2 = D_Tx / (2√2)
#
# Перевод в FWHM (формула 9):
#   θ_FWHM = √(ln2/2) · θ_e-2
#
# Усиление передающей антенны (формула 10):
#   g_Tx = 10·log10( (4·√2 / θ_e-2)² )
#        = 10·log10( (4·√(ln2) / θ_FWHM)² )
 
# Эффективный диаметр пучка с учётом коэффициента использования апертуры
space_D_e2_m = space_telescope_diam_m / math.sqrt(2)   # диаметр пучка на уровне 1/e²
space_omega0_m = space_D_e2_m / 2                       # радиус пучка ω0
 
theta_e2_rad   = 2 * wavelength_m / (math.pi * space_omega0_m)           # полный угол расходимости 1/e²
theta_FWHM_rad = math.sqrt(math.log(2) / 2) * theta_e2_rad               # угол расходимости FWHM
 
# Усиление передающей антенны (Giggenbach, формула 10)
g_Tx_dB = 10 * math.log10((4 * math.sqrt(2) / theta_e2_rad) ** 2)
 
# Внутренние потери передатчика (КПД оптического тракта)
a_Tx_dB = 10 * math.log10(space_efficiency)             # отрицательное значение
 
print(f"--- Передающий телескоп КА ---")
print(f"диаметр телескопа: {space_telescope_diam_m*1e3:.1f} мм")
print(f"D_e-2 (диаметр пучка 1/e²): {space_D_e2_m*1e3:.2f} мм")
print(f"ω0 (радиус пучка): {space_omega0_m*1e3:.2f} мм")
print(f"θ_e-2 (угол расходимости 1/e²): {theta_e2_rad*1e6:.2f} мкрад")
print(f"θ_FWHM (угол расходимости FWHM): {theta_FWHM_rad*1e6:.2f} мкрад")
print(f"усиление передающей антенны g_Tx: {g_Tx_dB:.3f} дБ")
print(f"внутренние потери передатчика a_Tx: {a_Tx_dB:.3f} дБ")
print()
 
# ── 5. Потери на ошибку наведения — формула (11)–(12) Giggenbach ─────────────
#
# Параметр β (соотношение угла расходимости и ошибки наведения):
#   β = (θ_FWHM)² / (4·ln2·σ_BW²)
#
# Средние потери на ошибку наведения:
#   a_BW = 10·log10( β / (β+1) )
 
beta_BW = (theta_FWHM_rad ** 2) / (4 * math.log(2) * space_pointing_error_sigma_rad ** 2)
a_BW_dB = 10 * math.log10(beta_BW / (beta_BW + 1))
 
print(f"--- Потери на ошибку наведения КА ---")
print(f"σ_BW (ошибка наведения КА): {space_pointing_error_sigma_rad*1e6:.2f} мкрад")
print(f"β (параметр бета-распределения): {beta_BW:.2f}")
print(f"a_BW (потери на ошибку наведения): {a_BW_dB:.3f} дБ")
print()
 
# ── 6. Приёмный телескоп наземной станции — формула (18) Giggenbach ──────────
#
#   g_Rx = 10·log10( 4π·A_Rx / λ² )
# где A_Rx — эффективная площадь апертуры приёмного телескопа.
#
# Для телескопа Кассегрена учитывается центральное экранирование
# через коэффициент использования апертуры gs_aperture_usage_coeff.
 
gs_geom_area_m2 = math.pi * (gs_telescope_diam_m ** 2) / 4
gs_eff_area_m2  = gs_aperture_usage_coeff * gs_geom_area_m2
 
g_Rx_dB = 10 * math.log10((4 * math.pi * gs_eff_area_m2) / (wavelength_m ** 2))
 
# Внутренние потери приёмного тракта (КПД телескопа + расщепитель и пр.)
a_Rx_dB = 10 * math.log10(gs_efficiency) + a_rx_internal_dB
 
print(f"--- Приёмный телескоп наземной станции ---")
print(f"диаметр телескопа НС: {gs_telescope_diam_m*1e2:.0f} см")
print(f"геометрическая площадь: {gs_geom_area_m2:.4f} м²")
print(f"эффективная площадь A_Rx: {gs_eff_area_m2:.4f} м²")
print(f"усиление приёмной антенны g_Rx: {g_Rx_dB:.3f} дБ")
print(f"внутренние потери приёмного тракта a_Rx: {a_Rx_dB:.3f} дБ")
print()
 
# ── 7. Суммарный бюджет линии — формула (1) Giggenbach ───────────────────────
#
#   p_Rx = p_Tx + a_Tx + g_Tx + a_BW + a_FSL + a_Atm + a_Sci + g_Rx + a_Rx
#
# Все величины в дБ/дБм; усиления положительны, потери отрицательны.
 
p_Rx_dBm = (space_laser_tx_power_dBm
            + a_Tx_dB
            + g_Tx_dB
            + a_BW_dB
            + a_FSL_dB
            + a_Atm_dB
            + a_sci_dB
            + g_Rx_dB
            + a_Rx_dB)
 
print(f"--- Суммарный бюджет ---")
print(f"мощность передатчика p_Tx:          {space_laser_tx_power_dBm:+.3f} дБм")
print(f"внутренние потери передатчика a_Tx:  {a_Tx_dB:+.3f} дБ")
print(f"усиление передающей антенны g_Tx:   {g_Tx_dB:+.3f} дБ")
print(f"потери на ошибку наведения a_BW:    {a_BW_dB:+.3f} дБ")
print(f"потери в своб. пространстве a_FSL:  {a_FSL_dB:+.3f} дБ")
print(f"атмосферное ослабление a_Atm:       {a_Atm_dB:+.3f} дБ")
print(f"потери на замирания a_Sci:           {a_sci_dB:+.3f} дБ")
print(f"усиление приёмной антенны g_Rx:     {g_Rx_dB:+.3f} дБ")
print(f"внутренние потери приёмника a_Rx:   {a_Rx_dB:+.3f} дБ")
print(f"мощность на детекторе p_Rx:         {p_Rx_dBm:+.3f} дБм")
print()
 
# ── 8. Требуемая мощность сигнала — формула (20) Giggenbach ──────────────────
#
# Для оптимизированного InGaAs-APD при BER=1e-3 (достаточно для FEC):
#   P_1E-3 = N_ppb · R · h·c/λ
# где:
#   N_ppb — число фотонов на бит
#   R     — скорость передачи данных, бит/с
#   h     — постоянная Планка
#   c     — скорость света
#   λ     — длина волны
h_planck = 6.62607015e-34           # постоянная Планка, Дж·с

data_rate_bps = data_flow_megabit_s * 1e6           # скорость передачи данных, бит/с
energy_per_photon_J = h_planck * speed_of_light_m_s / wavelength_m   # энергия фотона, Дж
 
P_required_W   = photons_per_bit * data_rate_bps * energy_per_photon_J
P_required_mW  = P_required_W * 1e3
P_required_dBm = 10 * math.log10(P_required_mW)
 
print(f"--- Требуемая мощность сигнала ---")
print(f"скорость передачи данных: {data_rate_bps/1e6:.3f} Мбит/с")
print(f"энергия фотона при λ=1550 нм: {energy_per_photon_J:.4e} Дж")
print(f"число фотонов на бит (N_ppb): {photons_per_bit} (для BER=1e-3, InGaAs-APD)")
print(f"требуемая мощность P_1E-3: {P_required_W*1e9:.4f} нВт  ({P_required_dBm:.3f} дБм)")
print()
 
# ── 9. Энергетический запас линии ────────────────────────────────────────────
 
margin_dB = p_Rx_dBm - P_required_dBm
 
print(f"--- Энергетический запас ---")
print(f"мощность на детекторе:   {p_Rx_dBm:.3f} дБм")
print(f"требуемая мощность:      {P_required_dBm:.3f} дБм")
print(f"запас энергетики линии:  {margin_dB:.3f} дБ")
print()
 
if margin_dB > 0:
    print("✅ Запас положительный — связь возможна")
else:
    print("❌ Запас отрицательный — требуется увеличение мощности,")
    print("   увеличение апертуры телескопа или уменьшение потерь")