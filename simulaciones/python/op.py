import numpy as np
import matplotlib.pyplot as plt

# --- Parámetros físicos de la planta ---
g = 9.81
a1 = 6.1360e-04
a2 = 5.4678e-04        
a3 = 3.5752e-04
a4 = 4.1820e-04 
k = 0.0017
gamma1 = 0.455
gamma2 = 0.4195

H_MAX = 0.60 # Altura máxima física de los estanques (60 cm)

# Matriz de Ganancia Estática M
M = np.array([
    [gamma1 * k, (1 - gamma2) * k],
    [(1 - gamma1) * k, gamma2 * k]
])
M_inv = np.linalg.inv(M)

# --- LÍMITES DE FLUJO ---
U_MIN = 0.26  # 26% Zona Muerta (Fricción estática)
U_MAX = 1.00  # 100% Saturación máxima

def calc_H(F, area):
    """Calcula la altura dado un caudal y el área de descarga"""
    return np.where(F < 0, -1, (F**2) / (2 * g * area**2))

def generar_poligono_principal(u_min, u_max, res=100):
    """Genera las coordenadas del polígono (Ambas bombas ON)"""
    u1_b = np.concatenate([np.linspace(u_min, u_max, res), np.full(res, u_max), np.linspace(u_max, u_min, res), np.full(res, u_min)])
    u2_b = np.concatenate([np.full(res, u_min), np.linspace(u_min, u_max, res), np.full(res, u_max), np.linspace(u_max, u_min, res)])
    
    F_b = M @ np.vstack([u1_b, u2_b])
    h1_b = calc_H(F_b[0, :], a1)
    h2_b = calc_H(F_b[1, :], a2)
    return h1_b, h2_b

def generar_trayectorias_1_bomba(u_min, u_max, res=100):
    """Genera las curvas utilizables cuando 1 sola bomba está encendida"""
    u_activa = np.linspace(u_min, u_max, res)
    
    # Cola 1: Bomba 1 OFF (0%), Bomba 2 ON
    F_cola1 = M @ np.vstack([np.zeros(res), u_activa])
    h1_c1 = calc_H(F_cola1[0, :], a1)
    h2_c1 = calc_H(F_cola1[1, :], a2)
    
    # Cola 2: Bomba 2 OFF (0%), Bomba 1 ON
    F_cola2 = M @ np.vstack([u_activa, np.zeros(res)])
    h1_c2 = calc_H(F_cola2[0, :], a1)
    h2_c2 = calc_H(F_cola2[1, :], a2)
    
    return (h1_c1, h2_c1), (h1_c2, h2_c2)

def es_alcanzable(u, tolerancia=1e-4):
    """Valida si un esfuerzo u es físicamente posible (0 o entre U_MIN y U_MAX)"""
    if abs(u) < tolerancia:
        return True
    if U_MIN - tolerancia <= u <= U_MAX + tolerancia:
        return True
    return False

def evaluar_puntos_utilizables(h1_sp, h2_sp):
    # 1. Calcular flujos requeridos
    F1_sp = a1 * np.sqrt(2 * g * h1_sp)
    F2_sp = a2 * np.sqrt(2 * g * h2_sp)
    
    # 2. Calcular esfuerzo requerido de bombas
    U_req = M_inv @ np.array([F1_sp, F2_sp])
    u1, u2 = U_req[0], U_req[1]
    
    # Limpieza numérica de ruidos cercanos a cero
    if abs(u1) < 1e-4: u1 = 0.0
    if abs(u2) < 1e-4: u2 = 0.0
    
    # 3. Validar si el punto es utilizable físicamente
    is_feasible = es_alcanzable(u1) and es_alcanzable(u2)
    
    # Evaluar rebalses superiores
    u1_eff = np.clip(u1, 0, 1)
    u2_eff = np.clip(u2, 0, 1)
    h3_calc = calc_H((1 - gamma2) * k * u2_eff, a3)
    h4_calc = calc_H((1 - gamma1) * k * u1_eff, a4)
    rebalse = h3_calc > H_MAX or h4_calc > H_MAX

    # --- GRAFICACIÓN DEL ESPACIO UTILIZABLE ---
    h1_poly, h2_poly = generar_poligono_principal(U_MIN, U_MAX)
    (h1_c1, h2_c1), (h1_c2, h2_c2) = generar_trayectorias_1_bomba(U_MIN, U_MAX)
    
    fig, ax = plt.subplots(figsize=(9, 7))
    
    # Rellenar Zona Principal (Ambas bombas)
    ax.fill(h1_poly, h2_poly, color='royalblue', alpha=0.3, label=f'Operación MIMO (Ambas Bombas {U_MIN*100:.0f}%-{U_MAX*100:.0f}%)')
    ax.plot(h1_poly, h2_poly, 'b-', linewidth=1.5)
    
    # Dibujar trayectorias degradadas (1 Bomba)
    ax.plot(h1_c1, h2_c1, color='purple', linestyle='--', linewidth=2.5, label='Operación SIMO (Bomba 1 OFF)')
    ax.plot(h1_c2, h2_c2, color='green', linestyle='--', linewidth=2.5, label='Operación SIMO (Bomba 2 OFF)')
    
    # Marcar el Origen
    ax.plot(0, 0, 'ko', markersize=8, label='Planta Apagada')

    # Evaluar color del SP solicitado
    if not is_feasible:
        color_sp, texto_estado = 'red', "INALCANZABLE (Saturación/Zona Muerta)"
    elif rebalse:
        color_sp, texto_estado = 'orange', "PELIGRO (Rebalse Superior)"
    elif (u1 == 0 and u2 != 0) or (u2 == 0 and u1 != 0):
        color_sp, texto_estado = 'purple', "ALCANZABLE (Operando con 1 Bomba)"
    else:
        color_sp, texto_estado = 'green', "ALCANZABLE (Operación MIMO Segura)"

    ax.plot(h1_sp, h2_sp, marker='X', markersize=12, color=color_sp, zorder=5, markeredgecolor='black')
    
    # Configuraciones visuales
    ax.set_title(f'Espacio de Puntos Utilizables\nSP = ({h1_sp*100:.1f} cm, {h2_sp*100:.1f} cm) | {texto_estado}', fontweight='bold')
    ax.set_xlabel('Nivel Estanque 1 [metros]')
    ax.set_ylabel('Nivel Estanque 2 [metros]')
    ax.set_xlim(-0.02, 0.6)
    ax.set_ylim(-0.02, 0.6)
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.7)
    
    plt.show()

    # REPORTE
    print("="*50)
    print(f"Setpoints Solicitados : h1 = {h1_sp*100:.1f} cm  |  h2 = {h2_sp*100:.1f} cm")
    print(f"Estado del Sistema    : {texto_estado}")
    print("-" * 50)
    print("ESFUERZOS REQUERIDOS (BOMBAS):")
    print(f"  > Bomba 1 (U1)      : {u1*100:.1f} %")
    print(f"  > Bomba 2 (U2)      : {u2*100:.1f} %")

# PRUEBA: Intenta usar un SP que caiga justo en la trayectoria de 1 bomba
# Por ejemplo, pide flujos proporcionales a las matrices para u1=0, u2=50%
u2_test = 0.50
F1_test = M[0,1]*u2_test
F2_test = M[1,1]*u2_test
h1_test = calc_H(F1_test, a1)
h2_test = calc_H(F2_test, a2)

print("Ingresa tu Setpoint a evaluar:")
sp_1 = float(input("Nivel L1 (metros) [ej. 0.013]: "))
sp_2 = float(input("Nivel L2 (metros) [ej. 0.054]: "))
evaluar_puntos_utilizables(sp_1*0.6/100, sp_2*0.6/100)