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


#=============================================================================
# II. БЮДЖЕТ ЛАЗЕРНОЙ ЛИНИИ СВЯЗИ КА–ЗЕМЛЯ
#
# [G] Giggenbach et al. — Link Budget Calculation in Optical LEO Satellite
#     Downlinks with On/Off-Keying and Large Signal Divergence, 2023
# [L] Liang et al. — Link Budget Analysis for Free-Space Optical Satellite
#     Networks, IEEE WoWMoM 2022
# =============================================================================
print()
print("=" * 60)
print("II. БЮДЖЕТ ЛАЗЕРНОЙ ЛИНИИ СВЯЗИ КА–ЗЕМЛЯ")
print("=" * 60)
 
# ── Исходные параметры ────────────────────────────────────────────────────────
epsilon_min_deg = 5.0
epsilon_min_rad = math.radians(epsilon_min_deg)
 
space_laser_tx_power_W   = 1.0
space_laser_tx_power_mW  = space_laser_tx_power_W * 1e3
space_laser_tx_power_dBm = 10 * math.log10(space_laser_tx_power_mW)
 
space_telescope_diam_m   = 5e-2       # диаметр передающего телескопа КА, м
gs_telescope_diam_m      = 0.5        # диаметр приёмного телескопа НС, м
 
space_efficiency           = 0.85
gs_efficiency              = 0.80
space_aperture_usage_coeff = 0.60
gs_aperture_usage_coeff    = 0.60
 
space_pointing_sigma_rad = 5e-6       # σ ошибки наведения КА, рад
 
a_rx_internal_dB = -4.1              # потери внутри НС (расщепитель и т.д.), дБ [G, Табл. V]
photons_per_bit  = 250               # фотонов/бит для InGaAs-APD, BER=1e-3 [G, ур. 20]
 
# ── Геометрия: максимальная дальность [G, ур. 14] ────────────────────────────
RE = r_earth_km
H0 = altitude_km
L_max_km = (math.sqrt((RE * math.sin(epsilon_min_rad)) ** 2 + 2 * H0 * RE + H0 ** 2)
            - RE * math.sin(epsilon_min_rad))
L_max_m  = L_max_km * 1e3
 
print(f"\n  Мощность передатчика КА: {space_laser_tx_power_dBm:.1f} дБм ({space_laser_tx_power_W:.0f} Вт)")
print(f"  Минимальный угол места:  {epsilon_min_deg}°")
print(f"  Максимальная дальность:  {L_max_km:.3f} км")
 
# ── Потери в свободном пространстве [G, ур. 13]: a_FSL = 20·log10(λ/(4πL)) ──
a_FSL_dB = 20 * math.log10(wavelength_m / (4 * math.pi * L_max_m))
 
# ── Передающая антенна КА [G, ур. 8–10] ─────────────────────────────────────
# Соотношение D_Tx = √2·D_e-2 (раздел III.A)
space_D_e2_m   = space_telescope_diam_m / math.sqrt(2)
space_omega0_m = space_D_e2_m / 2
theta_e2_rad   = 2 * wavelength_m / (math.pi * space_omega0_m)       # [G, ур. 8]
theta_FWHM_rad = math.sqrt(math.log(2) / 2) * theta_e2_rad           # [G, ур. 9]
g_Tx_dB        = 10 * math.log10((4 * math.sqrt(2) / theta_e2_rad) ** 2)  # [G, ур. 10]
a_Tx_dB        = 10 * math.log10(space_efficiency)
 
print(f"\n  Передающий телескоп КА (D={space_telescope_diam_m*1e3:.0f} мм):")
print(f"    θ_e-2 = {theta_e2_rad*1e6:.2f} мкрад,  θ_FWHM = {theta_FWHM_rad*1e6:.2f} мкрад")
print(f"    g_Tx  = {g_Tx_dB:.3f} дБ,  a_Tx = {a_Tx_dB:.3f} дБ")
 
# ── Потери на ошибку наведения [G, ур. 11–12] ────────────────────────────────
beta_BW = (theta_FWHM_rad ** 2) / (4 * math.log(2) * space_pointing_sigma_rad ** 2)
a_BW_dB = 10 * math.log10(beta_BW / (beta_BW + 1))
print(f"    β = {beta_BW:.2f},  a_BW = {a_BW_dB:.3f} дБ")
 
# ── Приёмная антенна НС [G, ур. 18]: g_Rx = 10·log10(4π·A_Rx/λ²) ───────────
gs_geom_area_m2 = math.pi * (gs_telescope_diam_m ** 2) / 4
gs_eff_area_m2  = gs_aperture_usage_coeff * gs_geom_area_m2
g_Rx_dB = 10 * math.log10((4 * math.pi * gs_eff_area_m2) / (wavelength_m ** 2))
a_Rx_dB = 10 * math.log10(gs_efficiency) + a_rx_internal_dB
print(f"\n  Приёмный телескоп НС (D={gs_telescope_diam_m*100:.0f} см):")
print(f"    A_eff = {gs_eff_area_m2:.4f} м²,  g_Rx = {g_Rx_dB:.3f} дБ,  a_Rx = {a_Rx_dB:.3f} дБ")
 
# ── Требуемая мощность приёмника [G, ур. 20]: P = N_ppb·R·hc/λ ─────────────
h_planck = 6.62607015e-34  # постоянная Планка, Дж·с
data_rate_bps     = data_flow_megabit_s * 1e6
energy_per_photon = h_planck * speed_of_light_m_s / wavelength_m
P_required_W      = photons_per_bit * data_rate_bps * energy_per_photon
P_required_mW     = P_required_W * 1e3
P_required_dBm    = 10 * math.log10(P_required_mW)
print(f"\n  Скорость передачи данных: {data_flow_megabit_s:.3f} Мбит/с")
print(f"  Требуемая мощность (250 фот/бит, BER=1e-3): {P_required_W*1e9:.3f} нВт  ({P_required_dBm:.3f} дБм)")
print(f"  Потери в своб. пространстве a_FSL: {a_FSL_dB:.3f} дБ")
 
# =============================================================================
# III. СРАВНЕНИЕ МОДЕЛЕЙ АТМОСФЕРНОГО ОСЛАБЛЕНИЯ
# =============================================================================
print()
print("=" * 60)
print("III. СРАВНЕНИЕ МОДЕЛЕЙ АТМОСФЕРНОГО ОСЛАБЛЕНИЯ")
print("=" * 60)
print(f"  Угол места: {epsilon_min_deg}° (наихудший случай)")
 
# ── Модель A: Giggenbach [G, ур. 16] ─────────────────────────────────────────
# Основана на зенитном коэффициенте пропускания Tz.
# Формула: a_Atm = 10·log10( Tz^(1/sin(ε)) )
# Tz=0.891 — наихудший сценарий из Таблицы II: λ=1550 нм, уровень моря,
# тропическая/городская среда. Учитывает только аэрозольное рассеяние
# и молекулярное поглощение (без облаков и тумана).
Tz = 0.891
a_Atm_G_dB = 10 * math.log10(Tz ** (1.0 / math.sin(epsilon_min_rad)))
 
print(f"\n  --- Модель A: Giggenbach [G, ур.16] ---")
print(f"  Параметр: Tz = {Tz}  (уровень моря, тропики, λ=1550нм, Таблица II)")
print(f"  Учитывает: аэрозольное и молекулярное ослабление (ясная погода)")
print(f"  a_Atm (Giggenbach) при ε={epsilon_min_deg}°:   {a_Atm_G_dB:.4f} дБ")
 
# ── Модель B: Liang et al. [L, ур. 6–11] ─────────────────────────────────────
# Учитывает два механизма рассеяния:
#   Im — рассеяние Ми (крупные частицы воды, нижняя атмосфера)
#   Ig — геометрическое рассеяние (туман, плотные облака)
# Суммарное ослабление: LA = Im · Ig  [L, ур. 11]
 

wavelength_nm = wavelength_m * 1e9        # λ в нм для эмпирических коэффициентов [L, ур. 6]
lam = wavelength_m * 1e6          # λ в мкм для эмпирических коэффициентов [L, ур. 6]
hE  = 0.0                    # высота НС над уровнем моря, км
hA  = 20.0                   # высота тропосферы, км
LW  = 0.5                    # содержание жидкой воды в облаке, г/м³
N   = 3.128e-4               # концентрация облачных частиц, см⁻³ [L, Табл. 3]
phi = 1.6                    # показатель размера частиц (модель Кима) [L, ур. 9]
 
# Эмпирические коэффициенты рассеяния Ми [L, ур. 6]:
a_c =  -0.000545 * lam**2 + 0.002  * lam - 0.0038
b_c =   0.00628  * lam**2 - 0.0232 * lam + 0.00439
c_c =  -0.028    * lam**2 + 0.101  * lam - 0.18
d_c =  -0.228    * lam**3 + 0.922  * lam**2 - 1.26 * lam + 0.719
 
rho  = a_c * hE**3 + b_c * hE**2 + c_c * hE + d_c   # коэффициент экстинкции
 
# Рассеяние Ми [L, ур. 7]: Im = exp(−ρ / sin(θE))
Im    = math.exp(-rho / math.sin(epsilon_min_rad))
Im_dB = 10 * math.log10(Im)
 
# Геометрическое рассеяние [L, ур. 8–10]:
V       = 1.002 / (LW * N) ** 0.6473              # видимость, км [L, ур. 8]
theta_A = (3.91 / V) * (wavelength_nm / 550) ** (-phi)  # коэф. ослабления, км⁻¹ [L, ур. 9]
dA      = (hA - hE) / math.sin(epsilon_min_rad)   # длина пути через тропосферу, км
Ig      = math.exp(-theta_A * dA)                 # [L, ур. 10]
Ig_dB   = 10 * math.log10(Ig)
 
# Суммарное атмосферное ослабление [L, ур. 11]:
LA          = Im * Ig
a_Atm_L_dB = 10 * math.log10(LA)
 
print(f"\n  --- Модель B: Liang et al. [L, ур.6–11] ---")
print(f"  Параметры: hE={hE} км, hA={hA} км, LW={LW} г/м³, N={N:.3e} см⁻³, φ={phi}")
print(f"  Учитывает: рассеяние Ми + геометрическое рассеяние (облака/туман)")
print(f"  ρ (коэф. экстинкции Ми):           {rho:.6f}")
print(f"  V (видимость):                     {V:.2f} км")
print(f"  θA (коэф. геом. рассеяния):        {theta_A:.6e} км⁻¹")
print(f"  dA (длина пути через атмосферу):   {dA:.2f} км")
print(f"  Im (потери рассеяния Ми):          {Im_dB:.4f} дБ")
print(f"  Ig (потери геом. рассеяния):       {Ig_dB:.4f} дБ")
print(f"  a_Atm (Liang) при ε={epsilon_min_deg}°:      {a_Atm_L_dB:.4f} дБ")
 
delta_dB = a_Atm_L_dB - a_Atm_G_dB
print(f"\n  ┌──────────────────────────────────────────────────────┐")
print(f"  │ СРАВНЕНИЕ МОДЕЛЕЙ при ε_min = {epsilon_min_deg}°                   │")
print(f"  │  Giggenbach (ясная погода, аэрозоль):  {a_Atm_G_dB:+6.2f} дБ  │")
print(f"  │  Liang (Ми + облачность/туман):        {a_Atm_L_dB:+6.2f} дБ  │")
print(f"  │  Разница (Liang − Giggenbach):         {delta_dB:+6.2f} дБ  │")
print(f"  │                                                      │")
print(f"  │  Модель Liang консервативнее на {abs(delta_dB):.2f} дБ:          │")
print(f"  │  дополнительно учитывает рассеяние Ми на водяных     │")
print(f"  │  каплях и геометрическое рассеяние в облаках/тумане. │")
print(f"  │  Модель Giggenbach — для условий ясного неба.        │")
print(f"  └──────────────────────────────────────────────────────┘")
 
# =============================================================================
# IV. ИТОГОВЫЙ БЮДЖЕТ ЛИНИИ — оба варианта ослабления
# =============================================================================
print()
print("=" * 60)
print("IV. ИТОГОВЫЙ БЮДЖЕТ ЛИНИИ СВЯЗИ [G, ур. 1]")
print("    p_Rx = p_Tx + a_Tx + g_Tx + a_BW + a_FSL + a_Atm + a_Sci + g_Rx + a_Rx")
print("=" * 60)
 
a_sci_dB = 0.0  # Сцинтилляция не учитываются при усреднении 100 мс [G, раздел VI]
 
def print_budget(label, a_Atm_dB):
    p_Rx_dBm  = (space_laser_tx_power_dBm + a_Tx_dB + g_Tx_dB + a_BW_dB
                 + a_FSL_dB + a_Atm_dB + a_sci_dB + g_Rx_dB + a_Rx_dB)
    margin_dB = p_Rx_dBm - P_required_dBm
    sym        = "✅" if margin_dB > 0 else "❌"
    print(f"\n  ── Модель атмосферы: {label} ──")
    print(f"  {'p_Tx  мощность передатчика':<40} {space_laser_tx_power_dBm:+8.3f} дБм")
    print(f"  {'a_Tx  потери в тракте КА':<40} {a_Tx_dB:+8.3f} дБ")
    print(f"  {'g_Tx  усиление антенны КА':<40} {g_Tx_dB:+8.3f} дБ")
    print(f"  {'a_BW  ошибка наведения (β-распред.)':<40} {a_BW_dB:+8.3f} дБ")
    print(f"  {'a_FSL потери в своб. пространстве':<40} {a_FSL_dB:+8.3f} дБ")
    print(f"  {'a_Atm атмосферное ослабление':<40} {a_Atm_dB:+8.3f} дБ")
    print(f"  {'a_Sci замирания (IRT)':<40} {a_sci_dB:+8.3f} дБ")
    print(f"  {'g_Rx  усиление антенны НС':<40} {g_Rx_dB:+8.3f} дБ")
    print(f"  {'a_Rx  потери в тракте НС':<40} {a_Rx_dB:+8.3f} дБ")
    print(f"  {'─'*50}")
    print(f"  {'p_Rx  мощность на детекторе':<40} {p_Rx_dBm:+8.3f} дБм")
    print(f"  {'P_req требуемая мощность':<40} {P_required_dBm:+8.3f} дБм")
    print(f"  {'ЗАПАС ЭНЕРГЕТИКИ ЛИНИИ':<40} {margin_dB:+8.3f} дБ  {sym}")
    return p_Rx_dBm, margin_dB
 
pRx_G, margin_G = print_budget("Giggenbach (ясная погода)", a_Atm_G_dB)
pRx_L, margin_L = print_budget("Liang (облачность/туман)", a_Atm_L_dB)
 
# =============================================================================
# V. АНАЛИЗ СЕАНСА ОПТИЧЕСКОЙ СВЯЗИ
# =============================================================================
print()
print("=" * 60)
print("V. АНАЛИЗ СЕАНСА ОПТИЧЕСКОЙ СВЯЗИ (DOWNLINK)")
print("=" * 60)
 
# ── 5.1 Геометрия и продолжительность сеанса ─────────────────────────────────
#
# КА совершает ровно 16 витков/сутки → трасса полностью повторяется ежедневно.
# Расстояние между соседними трассами: 360°/16 = 22.5° ≈ 2502 км (по экватору),
# что в ~85 раз больше ширины съёмочной полосы (29.5 км) — перекрытие отсутствует.
#
# Наземная станция принимает сигнал, пока КА выше ε_min = 5°.
# Геометрия из треугольника КА–НС–центр Земли [G, рис. 3]:
#   η_max   = arcsin( RE/(RE+H0) · cos(ε_min) )   — макс. надирный угол
#   λ_max   = π/2 − η_max − ε_min                  — макс. центральный угол
#   λ_max_km = λ_max · RE                           — дуга на поверхности Земли
#
# КА проходит дугу 2·λ_max_km вдоль трассы со скоростью v_shadow:
#   t_сеанс = 2 · λ_max_km / v_shadow
 
eta_max_rad    = math.asin((RE / (RE + H0)) * math.cos(epsilon_min_rad))
lambda_max_rad = math.pi / 2 - eta_max_rad - epsilon_min_rad
lambda_max_km  = lambda_max_rad * RE
t_session_s    = 2 * lambda_max_km / v_shadow_km_s
 
track_sep_km = math.radians(360.0 / cycles) * RE
 
print(f"\n  5.1 Геометрия орбиты")
print(f"  Высота орбиты:                    {H0:.2f} км")
print(f"  Период обращения:                 {period_s:.2f} с  ({period_s/60:.3f} мин)")
print(f"  Витков в сутки:                   {cycles}")
print(f"  Расстояние между трассами:        {track_sep_km:.1f} км")
print(f"  Ширина съёмочной полосы:          {fov_total_length_lambda_km:.2f} км")
print(f"  Перекрытие трасс:                 отсутствует ({track_sep_km/fov_total_length_lambda_km:.0f}× больше полосы)")
 
print(f"\n  5.2 Длительность сеанса связи")
print(f"  η_max (макс. надирный угол):      {math.degrees(eta_max_rad):.4f}°")
print(f"  λ_max (макс. центр. угол):        {math.degrees(lambda_max_rad):.4f}°  ({lambda_max_km:.2f} км)")
print(f"  Дуга видимости по Земле 2·λ_max:  {2*lambda_max_km:.2f} км")
print()
print(f"  ╔══════════════════════════════════════════════╗")
print(f"  ║  ДЛИТЕЛЬНОСТЬ СЕАНСА: {t_session_s:.1f} с ({t_session_s/60:.3f} мин)   ║")
print(f"  ╚══════════════════════════════════════════════╝")
print(f"  (Типовое значение для LEO по Giggenbach: «несколько–10 мин» ✓)")
 
# ── 5.2 Освещённость: доля орбиты над lit-стороной Земли ─────────────────────
#
# Для круговой орбиты при нулевом бета-угле (наихудший случай):
#   f_eclipse = arccos( sqrt(H0²+2·RE·H0) / (RE+H0) ) / π
#
# Это стандартная геометрическая формула для цилиндрической тени Земли.
# При β=0° тень максимальна; реальное значение меньше.
 
f_eclipse  = math.acos(math.sqrt(H0**2 + 2*RE*H0) / (RE + H0)) / math.pi
f_sunlit   = 1.0 - f_eclipse
t_sunlit_s = f_sunlit * period_s
t_eclipse_s = f_eclipse * period_s
 
print(f"\n  5.3 Освещённость (бета-угол = 0°, наихудший случай по тени)")
print(f"  Доля в тени:      {f_eclipse*100:.2f}%  ({t_eclipse_s:.1f} с = {t_eclipse_s/60:.2f} мин)")
print(f"  Доля на свету:    {f_sunlit*100:.2f}%  ({t_sunlit_s:.1f} с = {t_sunlit_s/60:.2f} мин)")
 
# ── 5.3 Накопление и передача данных ─────────────────────────────────────────
#
# КА непрерывно снимает, пока:
#   (а) находится над освещённой стороной Земли
#   (б) земля попадает в съёмочную полосу шириной fov_total_length_lambda_km
#
# Условие (б) выполняется всегда, пока работает камера (это ограничение
# на ширину полосы, а не на продолжительность съёмки). Поэтому
# время активной съёмки = t_sunlit_s за каждый виток.
#
# Объём данных, НАКОПЛЕННЫХ за освещённую часть одного витка:
data_per_lit_pass_Mb = data_flow_megabit_s * t_sunlit_s
 
# Длина отснятой полосы за освещённую часть витка:
strip_length_km = v_shadow_km_s * t_sunlit_s
 
# Объём данных, ПЕРЕДАННЫХ за один сеанс (канал работает на полной скорости):
data_per_session_Mb = data_flow_megabit_s * t_session_s
 
# Суточное накопление (16 витков):
data_per_day_Mb = data_per_lit_pass_Mb * cycles
 
# ── 5.4 Суточный баланс связи ─────────────────────────────────────────────────
# Для полярной орбиты с повторяющейся трассой каждый наземный пункт
# попадает под трассу не чаще 1 раза в сутки. Однако за счёт угла
# видимости (2·λ_max ≈ 2766 км проекция на Земле) станция может
# принимать КА 2–4 раза в сутки для разных витков.
 
print(f"\n  5.4 Данные: накопление и передача")
print(f"  Поток ДЗЗ:                        {data_flow_megabit_s:.3f} Мбит/с")
print(f"  Длина снятой полосы (1 виток):    {strip_length_km:.1f} км  ×  {fov_total_length_lambda_km:.2f} км (ширина)")
print()
print(f"  За один виток:")
print(f"    — накоплено данных (свет. дуга):{data_per_lit_pass_Mb:>12.2f} Мбит  ({data_per_lit_pass_Mb/8e3:.3f} ГБ)")
print()
print(f"  ╔══════════════════════════════════════════════════════════╗")
print(f"  ║  ОДИН СЕАНС СВЯЗИ  (t_сеанс = {t_session_s:.1f} с = {t_session_s/60:.3f} мин):      ║")
print(f"  ║  Объём переданных данных:                                ║")
print(f"  ║    {data_per_session_Mb:>12.2f} Мбит                                   ║")
print(f"  ║  = {data_per_session_Mb/1e3:>12.3f} Гбит                                   ║")
print(f"  ║  = {data_per_session_Mb/8e3:>12.4f} ГБ                                     ║")
print(f"  ╚══════════════════════════════════════════════════════════╝")
 
ratio_pass_to_session = data_per_lit_pass_Mb / data_per_session_Mb
print(f"\n  Соотношение накоплено/сеанс:      {ratio_pass_to_session:.2f}×")
print(f"  → Один сеанс передаёт {1/ratio_pass_to_session*100:.1f}% данных одного освещённого витка")
print(f"  → Для полной разгрузки одного витка нужно {math.ceil(ratio_pass_to_session)} сеанса")
 
print(f"\n  За сутки ({cycles} витков):")
print(f"    — накоплено:                   {data_per_day_Mb/1e3:>10.3f} Гбит  ({data_per_day_Mb/8e3:.2f} ГБ)")
print()
print(f"  Суточный баланс (сценарии по числу сеансов/сут):")
print(f"  {'Сеансов/сут':<14} {'Передано, ГБ':<16} {'% от суточн. накопления'}")
print(f"  {'─'*55}")
for n_pass in [1, 2, 3, 4, 5]:
    tx_gb  = n_pass * data_per_session_Mb / 8e3
    pct    = tx_gb / (data_per_day_Mb / 8e3) * 100
    note   = " ← типовое" if n_pass == 3 else ""
    print(f"  {n_pass:<14} {tx_gb:<16.3f} {pct:.2f}%{note}")
 
# ── 5.5 Ресурс миссии ─────────────────────────────────────────────────────────
print(f"\n  5.5 Оценка объёма данных за ресурс миссии (6U CubeSat, 12 кг)")
print(f"  ({'3 сеанса/сут — типовое значение'})")
n_sessions_typ = 3
print(f"\n  {'Сценарий':<30} {'Накоплено, ГБ':>14} {'Передано, ГБ':>14}")
print(f"  {'─'*60}")
for label_m, days in [("Мин. ресурс на рабочей орб. (15 сут)", 15),
                      ("Макс. ресурс на рабочей орб. (30 сут)", 30),
                      ("Мин. до сход. с орбиты (65 сут)",       65),
                      ("Макс. до сход. с орбиты (130 сут)",    130)]:
    acc = days * data_per_day_Mb / 8e3
    tx  = days * n_sessions_typ * data_per_session_Mb / 8e3
    print(f"  {label_m:<38} {acc:>10.1f}     {tx:>10.1f}")
 
print()
print("  ВЫВОД: передача данных ограничена не пропускной способностью")
print(f"  лазерного канала ({data_flow_megabit_s:.0f} Мбит/с, запас по Giggenbach {margin_G:+.1f} дБ),")
print(f"  а числом сеансов связи в сутки и ёмкостью бортового накопителя.")
print(f"  Для полной передачи суточного объёма съёмки необходимо не менее")
print(f"  {math.ceil(data_per_day_Mb/data_per_session_Mb)} сеансов в сутки — что практически недостижимо")
print(f"  с одной наземной станцией. Рекомендуется сеть из нескольких НС.")