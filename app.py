import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, acos, radians, pi

def solar_declination(day_of_year):
    return 0.409 * sin((2 * pi / 365) * day_of_year - 1.39)

def sunset_hour_angle(lat_rad, delta):
    return acos(-tan(lat_rad) * tan(delta))

def extraterrestrial_radiation(lat_deg, day_of_year):
    lat_rad = radians(lat_deg)
    dr = 1 + 0.033 * cos((2 * pi / 365) * day_of_year)
    delta = solar_declination(day_of_year)
    ws = sunset_hour_angle(lat_rad, delta)
    Ra = (24 * 60 / pi) * 0.0820 * dr * (
        ws * sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * sin(ws)
    )
    return Ra

def eto_penman_monteith(T, RH, u2, Rs, P, albedo):
    delta = 4098 * (0.6108 * np.exp((17.27 * T) / (T + 237.3))) / ((T + 237.3) ** 2)
    gamma = 0.665e-3 * P
    es = 0.6108 * np.exp((17.27 * T) / (T + 237.3))
    ea = es * (RH / 100)
    Rns = (1 - albedo) * Rs
    Rnl = 4.903e-9 * (((T + 273.16)**4 + (T + 273.16)**4)/2) *           (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * min(1, Rs / (0.75 + 2e-5 * 101.3)) - 0.35)
    Rn = Rns - Rnl
    G = 0
    eto = (0.408 * delta * (Rn - G) + gamma * (900 / (T + 273)) * u2 * (es - ea)) /           (delta + gamma * (1 + 0.34 * u2))
    return max(0, eto)

st.set_page_config(layout="wide")
st.title("🌤️ Simulador de Evapotranspiración (ETo) - Penman-Monteith FAO 56")
st.markdown("Visualiza la variación de la ETo durante el año ajustando variables climáticas realistas.")

st.latex(r"ETo = rac{0.408 \cdot \Delta (R_n - G) + \gamma \cdot rac{900}{T + 273} \cdot u_2 (e_s - e_a)}{\Delta + \gamma (1 + 0.34 u_2)}")

st.markdown("""
**Unidades**:
- ETo: mm/día  
- T: °C  
- RH: %  
- u₂: m/s  
- P: kPa  
- α (albedo): adimensional  
""")

col1, col2 = st.columns([1, 2])

with col1:
    T = st.slider("Temperatura media (°C)", -10.0, 50.0, 25.0)
    RH = st.slider("Humedad relativa (%)", 5.0, 100.0, 60.0)
    u2 = st.slider("Velocidad del viento (m/s)", 0.1, 10.0, 2.0)
    P = st.slider("Presión atmosférica (kPa)", 60.0, 110.0, 101.3)

    cultivos = {
        "Pasto corto (referencia)": 0.23,
        "Maíz": 0.20,
        "Arroz": 0.25,
        "Caña de azúcar": 0.18,
        "Trigo": 0.23,
        "Superficie clara (arena/nieve)": 0.30,
        "Superficie oscura (suelo húmedo)": 0.10
    }
    cultivo = st.selectbox("Superficie/cultivo", list(cultivos.keys()))
    albedo = cultivos[cultivo]

    lat = st.number_input("Latitud (°)", -66.0, 66.0, 10.0)

with col2:
    days = np.arange(1, 366)
    Ra_series = [extraterrestrial_radiation(lat, d) for d in days]
    Rs_series = [0.75 * Ra for Ra in Ra_series]
    eto_series = [eto_penman_monteith(T, RH, u2, Rs, P, albedo) for Rs in Rs_series]

    fig, ax = plt.subplots()
    ax.plot(days, eto_series, label="ETo diaria (mm)")
    ax.set_title("Variación de ETo durante el año")
    ax.set_xlabel("Día del año")
    ax.set_ylabel("ETo (mm/día)")
    ax.grid(True)
    ax.legend()
    st.pyplot(fig)

    st.success(f"🌱 Valor promedio anual estimado de ETo: {np.mean(eto_series):.2f} mm/día")
