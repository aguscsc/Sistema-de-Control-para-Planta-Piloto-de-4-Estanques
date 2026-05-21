import numpy as np
import matplotlib.pyplot as plt

# --- Parámetros físicos de la planta ---
g = 9.81
a1 = 6.1360e-04; a2 = 5.4678e-04        # Áreas de fuga inferiores
a3 = 3.5752e-04; a4 = 4.1820e-04 # Áreas de fuga superiores
k = 0.0017
gamma1 = 0.455; gamma2 = 0.4195

H_MAX = 0.60 # Altura máxima física de los estanques (60 cm)

# Matriz de Ganancia Estática M
M = np.array([
    [gamma1 * k, (1 - gamma2) * k],
    [(1 - gamma1) * k, gamma2 * k]
])
M_inv = np.linalg.inv(M)

# --- LÍMITES  DE FLUJO---
U_MIN = 0.26  # 40% Zona Muerta, 26% de flujo
U_MAX = 1.00  # 100% Saturación máxima

def obtener_limites_reales(F1_sp, F2_sp):
    I11, I12 = M_inv[0, 0], M_inv[0, 1]
    I21, I22 = M_inv[1, 0], M_inv[1, 1]

    # Límites para F2 dado F1
    bounds_F2 = [
        ('u1_min', (U_MIN - I11*F1_sp)/I12, I12),
        ('u1_max', (U_MAX - I11*F1_sp)/I12, I12),
        ('u2_min', (U_MIN - I21*F1_sp)/I22, I22),
        ('u2_max', (U_MAX - I21*F1_sp)/I22, I22)
    ]
    F2_min, F2_max = 0.0, float('inf')
    for name, val, coeff in bounds_F2:
        if name.endswith('_min'):
            if coeff > 0: F2_min = max(F2_min, val)
            else: F2_max = min(F2_max, val)
        elif name.endswith('_max'):
            if coeff > 0: F2_max = min(F2_max, val)
            else: F2_min = max(F2_min, val)

    h2_min = (F2_min**2)/(2*g*a2**2) if F2_max >= F2_min and F2_min >= 0 else 0
    h2_max = (F2_max**2)/(2*g*a2**2) if F2_max >= F2_min and F2_max > 0 else -1

    # Límites para F1 dado F2
    bounds_F1 = [
        ('u1_min', (U_MIN - I12*F2_sp)/I11, I11),
        ('u1_max', (U_MAX - I12*F2_sp)/I11, I11),
        ('u2_min', (U_MIN - I22*F2_sp)/I21, I21),
        ('u2_max', (U_MAX - I22*F2_sp)/I21, I21)
    ]
    F1_min, F1_max = 0.0, float('inf')
    for name, val, coeff in bounds_F1:
        if name.endswith('_min'):
            if coeff > 0: F1_min = max(F1_min, val)
            else: F1_max = min(F1_max, val)
        elif name.endswith('_max'):
            if coeff > 0: F1_max = min(F1_max, val)
            else: F1_min = max(F1_min, val)

    h1_min = (F1_min**2)/(2*g*a1**2) if F1_max >= F1_min and F1_min >= 0 else 0
    h1_max = (F1_max**2)/(2*g*a1**2) if F1_max >= F1_min and F1_max > 0 else -1

    return h1_min, h1_max, h2_min, h2_max

def calcular_y_graficar(h1_sp, h2_sp):
    # 1. Calcular flujos inferiores requeridos
    F1_sp = a1 * np.sqrt(2 * g * h1_sp)
    F2_sp = a2 * np.sqrt(2 * g * h2_sp)
    
    # Calcular esfuerzos de bombas (Inversa)
    U = M_inv @ np.array([F1_sp, F2_sp])
    u1, u2 = U[0], U[1]
    
    # Evaluar condiciones de seguridad de las bombas
    is_feasible = (U_MIN <= u1 <= U_MAX) and (U_MIN <= u2 <= U_MAX)
    
    # Calcular el nivel resultante en Estanques 3 y 4 (Los Esclavos)
    h3_result = (1 / (2 * g)) * (((1 - gamma2) * k * u2) / a3)**2
    h4_result = (1 / (2 * g)) * (((1 - gamma1) * k * u1) / a4)**2
    
    # Comprobar si hay rebalse físico en los estanques superiores
    rebalse = (h3_result > H_MAX) or (h4_result > H_MAX)

    # Obtener los límites dinámicos para los maestros
    h1_min, h1_max, h2_min, h2_max = obtener_limites_reales(F1_sp, F2_sp)

    # --- IMPRIMIR LA HOJA DE DATOS DE OPERACIÓN ---
    print("="*50)
    print("   HOJA DE PUNTOS DE OPERACIÓN - PLANTA JOHANSSON")
    print("="*50)
    print(f"Setpoints Solicitados : h1 = {h1_sp*100:.1f} cm  |  h2 = {h2_sp*100:.1f} cm")
    
    if is_feasible and not rebalse:
        print("Estado del Sistema    : [OK] - FÍSICAMENTE ALCANZABLE")
    elif not is_feasible:
        print("Estado del Sistema    : [ERROR] - BOMBAS SATURADAS O EN ZONA MUERTA")
    elif rebalse:
        print("Estado del Sistema    : [PELIGRO] - REBALSE EN ESTANQUES SUPERIORES")
        
    print("-" * 50)
    print("ESFUERZOS REQUERIDOS (BOMBAS):")
    print(f"  > Bomba 1 (U1)      : {u1*100:.1f} %")
    print(f"  > Bomba 2 (U2)      : {u2*100:.1f} %")
    
    print("-" * 50)
    print("CONSECUENCIA EN ESTANQUES SUPERIORES:")
    print(f"  > Estanque 3 (h3)   : {h3_result*100:.1f} cm " + ("(¡REBALSE!)" if h3_result > H_MAX else ""))
    print(f"  > Estanque 4 (h4)   : {h4_result*100:.1f} cm " + ("(¡REBALSE!)" if h4_result > H_MAX else ""))
    
    print("-" * 50)
    print("RANGOS DE OPERACIÓN DINÁMICOS:")
    if is_feasible:
        print(f"  > Si mantiene h2 en {h2_sp*100:.1f} cm, h1 puede variar entre: [{h1_min*100:.1f} - {h1_max*100:.1f} cm]")
        print(f"  > Si mantiene h1 en {h1_sp*100:.1f} cm, h2 puede variar entre: [{h2_min*100:.1f} - {h2_max*100:.1f} cm]")
    else:
        print("  > Rangos no disponibles debido a saturación.")
    print("="*50)

calcular_y_graficar(13*0.6/100, 21*0.6/100) # Ejemplo de punto de operación (h1=25.17 cm, h2=27.5 cm) con ajuste al 60% para evitar zona muerta
