import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos de la planta
g = 9.81; a = 0.001; k = 0.003
gamma1 = 0.455; gamma2 = 0.4195

M = np.array([
    [gamma1 * k, (1 - gamma2) * k],
    [(1 - gamma1) * k, gamma2 * k]
])
M_inv = np.linalg.inv(M)

# --- LÍMITES FÍSICOS REALES ---
U_MIN = 0.40  # 40% Zona Muerta
U_MAX = 1.00  # 100% Saturación máxima

def obtener_limites_reales(F1_sp, F2_sp):
    I11, I12 = M_inv[0, 0], M_inv[0, 1]
    I21, I22 = M_inv[1, 0], M_inv[1, 1]

    # Límites de F2 dado F1 usando U_MIN (0.4) en vez de 0.0
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

    h2_min = (F2_min**2)/(2*g*a**2) if F2_max >= F2_min and F2_min >= 0 else 0
    h2_max = (F2_max**2)/(2*g*a**2) if F2_max >= F2_min and F2_max > 0 else -1

    # Límites de F1 dado F2 usando U_MIN (0.4) en vez de 0.0
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

    h1_min = (F1_min**2)/(2*g*a**2) if F1_max >= F1_min and F1_min >= 0 else 0
    h1_max = (F1_max**2)/(2*g*a**2) if F1_max >= F1_min and F1_max > 0 else -1

    return h1_min, h1_max, h2_min, h2_max

def generar_perimetro(u_min, u_max, res=200):
    """Genera las coordenadas del polígono para graficar la zona"""
    u1_b = np.concatenate([np.linspace(u_min,u_max,res), np.full(res, u_max), np.linspace(u_max,u_min,res), np.full(res, u_min)])
    u2_b = np.concatenate([np.full(res, u_min), np.linspace(u_min,u_max,res), np.full(res, u_max), np.linspace(u_max,u_min,res)])
    F_b = M @ np.vstack([u1_b, u2_b])
    return F_b**2 / (2 * g * a**2)

def graficar_punto_real(h1_sp, h2_sp):
    # Calcular esfuerzos requeridos
    F1_sp = a * np.sqrt(2 * g * h1_sp)
    F2_sp = a * np.sqrt(2 * g * h2_sp)
    U = M_inv @ np.array([F1_sp, F2_sp])
    u1, u2 = U[0], U[1]
    
    # Evaluar condiciones
    is_feasible_real = (U_MIN <= u1 <= U_MAX) and (U_MIN <= u2 <= U_MAX)
    is_feasible_teorico = (0.0 <= u1 <= U_MAX) and (0.0 <= u2 <= U_MAX)

    h1_min, h1_max, h2_min, h2_max = obtener_limites_reales(F1_sp, F2_sp)

    #  Generar polígonos
    H_ideal = generar_perimetro(0.0, 1.0)
    H_real = generar_perimetro(U_MIN, U_MAX)

    # Crear el Gráfico
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Zona Teórica (0% - 100%)
    ax.plot(H_ideal[0,:], H_ideal[1,:], color='gray', linestyle='--', linewidth=1.5, alpha=0.8)
    ax.fill(H_ideal[0,:], H_ideal[1,:], color='lightgray', alpha=0.3, label='Zona Teórica Ideal (0-100%)')
    
    # Zona Real (40% - 100%)
    ax.plot(H_real[0,:], H_real[1,:], 'k-', linewidth=2)
    ax.fill(H_real[0,:], H_real[1,:], color='royalblue', alpha=0.6, label=f'Zona Controlable Real ({U_MIN*100:.0f}%-100%)')

    # Evaluar semáforo de seguridad
    if is_feasible_real:
        c = 'green'; st = "ALCANZABLE (ZONA SEGURA)"
    elif is_feasible_teorico:
        c = 'orange'; st = "PELIGRO (BOMBA EN ZONA MUERTA)"
    else:
        c = 'red'; st = "INALCANZABLE (SATURACIÓN > 100%)"

    # Marcar el punto
    ax.plot(h1_sp, h2_sp, marker='X', markersize=12, color=c, zorder=5)

    # Dibujar líneas de rango si es posible
    if h1_max != -1 and is_feasible_real:
        ax.plot([h1_min, h1_max], [h2_sp, h2_sp], color=c, linestyle='-', linewidth=3, label=f"Rango h1: [{h1_min*100:.1f}, {h1_max*100:.1f}] cm")
    if h2_max != -1 and is_feasible_real:
        ax.plot([h1_sp, h1_sp], [h2_min, h2_max], color=c, linestyle='-', linewidth=3, label=f"Rango h2: [{h2_min*100:.1f}, {h2_max*100:.1f}] cm")

    # Configuraciones visuales
    ax.set_title(f'Envolvente Operativa con Zona Muerta (Deadband)\nSP1={h1_sp*100:.0f}cm, SP2={h2_sp*100:.0f}cm | {st}\nU1Req: {u1*100:.1f}% | U2Req: {u2*100:.1f}%', fontsize=13, fontweight='bold', color=c)
    ax.set_xlabel('Nivel Estanque 1 [metros]', fontsize=12)
    ax.set_ylabel('Nivel Estanque 2 [metros]', fontsize=12)
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 0.6)
    ax.legend(loc='upper right')
    ax.grid(True, linestyle=':', alpha=0.7)
    
    plt.show()

    # Reporte en Consola
    print(f"=== REPORTE DINÁMICO ===")
    print(f"Esfuerzo Bomba 1: {u1*100:.1f}%")
    print(f"Esfuerzo Bomba 2: {u2*100:.1f}%")
    print(f"Estado: {st}")
    print("-" * 25)
    if is_feasible_real:
        print(f"Rango utilizable para Estanque 1: [{h1_min*100:.1f} cm a {h1_max*100:.1f} cm]")
        print(f"Rango utilizable para Estanque 2: [{h2_min*100:.1f} cm a {h2_max*100:.1f} cm]")
    else:
        print("El punto cae fuera de la zona segura. No se puede garantizar el control.")

sp_1 = float(input("sp 1:"))
sp_2 = float(input("sp 2:"))
graficar_punto_real(sp_1, sp_2)
