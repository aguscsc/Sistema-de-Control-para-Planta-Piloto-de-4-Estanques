"""
stress_horext.py
=================
Herramienta generica de stress-test para los 4 lazos HOREXT del L5K
(F1, F2, L3, L4). Reemplaza / generaliza a stress_l3.py.

Replica, scan a scan (Tm=0.1s), el comportamiento EXACTO del ST:
  - PRBS (LFSR de 7 bits, taps 6/7, semilla 0011100)
  - Controlador (manual con rate-limit, o HOREXT Estrategia 1 con la
    recursion geometrica p/beta, o PID velocity-form) con divisor de scans
  - Para L3/L4: cascada hacia el lazo de flujo interno (F2/F1) con su propio
    PID + rate-limiter, y planta de flujo (polo rapido) + estanque (polo lento)
  - Para F1/F2: planta directa bomba->flujo (polo rapido) + rate-limiter
  - estimador RLS-Fortescue generico para orden n y horizonte T cualquiera

Dos modos de uso:
  1) identificar(cfg, ...)   -> corre el PRBS en lazo abierto, reporta a0,b,...
  2) lazo_cerrado(cfg, ...)  -> activa HOREXT con un escalon de SP, reporta
                                overshoot/undershoot, saturacion, oscilacion.

Simplificaciones conocidas (no modeladas):
  - Feedforward (FF_L3/FF_4): es una suma estatica que depende de SP; se omite
    porque este modelo trabaja en espacio "porcentaje" abstracto, no en
    unidades fisicas SI. No afecta la identificacion (FF no se mueve durante
    el PRBS de nivel si SP se mantiene fijo).
  - Interaccion cruzada entre estanques (gamma1/gamma2, matriz I): cada lazo
    se simula desacoplado. Introduce optimismo leve respecto a la planta real
    cuando ambas bombas trabajan simultaneamente.

Uso rapido (linea de comandos):
    python stress_horext.py ident --loop L3
    python stress_horext.py ident --loop L3 --no-noise
    python stress_horext.py ident --loop L3 --filtro
    python stress_horext.py step  --loop L3 --sp-step -10
    python stress_horext.py step  --loop F2 --sp-step 15 --adaptive
    python stress_horext.py bateria      # corre todo, todos los lazos
"""
import argparse
import dataclasses
import numpy as np

Tm = 0.1  # tiempo base de la tarea periodica, igual que en el L5K


# =============================================================================
# PRBS (LFSR de 7 bits, igual al de PRBS_1/2/L3/L4 en el L5K)
# =============================================================================
class PRBS:
    def __init__(self, delta, amp, cv0, seed=(0, 0, 1, 1, 1, 0, 0)):
        self.l = list(seed)
        self.delta = delta
        self.act = 0
        self.pmax = cv0 + amp / 2.0
        self.pmin = cv0 - amp / 2.0
        self.out = cv0

    def step(self):
        self.act -= 1
        if self.act <= 0:
            self.act = int(self.delta / Tm)
            s = self.l[5] ^ self.l[6]
            self.l = [s] + self.l[:6]
            self.out = self.pmax if s else self.pmin
        return self.out


# =============================================================================
# Estimador RLS-Fortescue generico (orden n, horizonte T cualquiera)
# =============================================================================
class RLSHorext:
    """NOTA: p_ceiling es una salvaguarda que el ST del PLC NO tiene. El
    'windup' de covarianza (P crece sin limite cuando Sigma es muy chico
    respecto al ruido real del sensor) es un riesgo real, no solo de esta
    simulacion -- vale la pena agregar un chequeo de Traza_P similar en el
    estimador real (ademas del que ya existe para Denom_K <= 0)."""

    def __init__(self, orden, horizonte, sigma, p_ceiling=1e5):
        self.orden = orden
        self.horizonte = horizonte
        self.n_beta = orden + horizonte - 1
        self.npar = 2 * orden + horizonte - 1
        self.theta = np.zeros(self.npar)
        self.P = np.eye(self.npar) * 10000.0
        self.dy_hist = np.zeros(16)
        self.du_hist = np.zeros(16)
        self.sigma = sigma
        self.p_ceiling = p_ceiling
        self.lam_hist = []
        self.resets_por_windup = 0

    def update(self, dy, du_realizado):
        n, nb = self.orden, self.n_beta
        phi = np.zeros(self.npar)
        phi[:n] = -self.dy_hist[:n]
        phi[n:n + nb] = self.du_hist[:nb]

        est = phi @ self.theta
        e_pred = dy - est
        Pphi = self.P @ phi
        lam_dem = 1.0 + phi @ Pphi
        if lam_dem > 1e-6 and self.sigma > 1e-6:
            lam = 1.0 - (e_pred ** 2) / (lam_dem * self.sigma)
            lam = min(max(lam, 0.5), 1.0)
        else:
            lam = 1.0

        den = lam + phi @ Pphi
        if den > 1e-6:
            K = Pphi / den
            self.theta = self.theta + K * e_pred
            self.P = (self.P - np.outer(K, Pphi)) / lam

        if np.trace(self.P) > self.p_ceiling or not np.all(np.isfinite(self.P)):
            self.P = np.eye(self.npar) * 10000.0
            self.resets_por_windup += 1

        self.dy_hist[1:] = self.dy_hist[:-1]
        self.dy_hist[0] = dy
        self.du_hist[1:] = self.du_hist[:-1]
        self.du_hist[0] = du_realizado
        self.lam_hist.append(lam)

    @property
    def a(self):
        return -self.theta[:self.orden]

    @property
    def b(self):
        return self.theta[self.orden:self.orden + self.n_beta]

    @property
    def traza_p(self):
        return float(np.trace(self.P))


# =============================================================================
# Ley de control HOREXT (Estrategia 1) -- misma recursion geometrica del ST
# =============================================================================
def horext_du(a_coefs, b_coefs, horizonte, dy, y_now, sp, ganancia_min,
              rate, nb_deadband=0.0):
    p = float(np.sum(a_coefs))
    b = float(np.sum(b_coefs))

    p_pow, a0_usado, beta_run, ganancia = 1.0, 0.0, 0.0, 0.0
    for _ in range(horizonte):
        p_pow *= p
        a0_usado += p_pow
        beta_run = p * beta_run + b
        ganancia += beta_run

    ganancia_usada = ganancia
    if ganancia_usada >= 0.0:
        ganancia_usada = max(ganancia_usada, ganancia_min)
    else:
        ganancia_usada = min(ganancia_usada, -ganancia_min)

    h_libre = y_now + a0_usado * dy
    e_pred = sp - h_libre
    if abs(e_pred) < nb_deadband:
        e_pred = 0.0

    du = e_pred / ganancia_usada
    du = float(np.clip(du, -rate, rate))
    return du, ganancia, a0_usado


def pid_du(Kc, Ki, Ts, error, error_prev, rate):
    du = Kc * (1.0 + Ki * Ts / 2.0) * error - Kc * (1.0 - Ki * Ts / 2.0) * error_prev
    return float(np.clip(du, -rate, rate))


# =============================================================================
# Configuracion de cada lazo (valores por defecto = parametros del L5K)
# =============================================================================
@dataclasses.dataclass
class LoopConfig:
    nombre: str
    kind: str            # "directo" (F1/F2) o "cascada" (L3/L4)
    input_scans: int      # divisor de scans -> Ts_lento = Tm*input_scans
    dumax: float          # %/s rate limiter propio del lazo
    orden: int
    horizonte: int
    sigma: float
    delta: float          # periodo de simbolo PRBS [s]
    amp: float            # amplitud PRBS [%]
    ganancia_min: float
    cv0: float             # punto de operacion (bomba o SP_flujo) [%]
    pv0: float             # punto de operacion (flujo o nivel medido) [%]
    tau_planta: float       # polo fisico rapido bomba->flujo [s]
    noise_pv: float          # sigma ruido sensor propio (flujo o nivel)
    # solo cascada:
    inner: "LoopConfig" = None
    tau_tanque: float = None
    noise_inner: float = None


# Lazo de flujo interno F2 (usado directo y como esclavo de L3)
F2 = LoopConfig(
    nombre="F2", kind="directo", input_scans=10, dumax=1.0,
    orden=1, horizonte=3, sigma=20.0, delta=10.0, amp=10.0,
    ganancia_min=0.01, cv0=55.0, pv0=55.0, tau_planta=1.8, noise_pv=0.15,
)
F1 = LoopConfig(
    nombre="F1", kind="directo", input_scans=10, dumax=1.0,
    orden=1, horizonte=3, sigma=20.0, delta=10.0, amp=10.0,
    ganancia_min=0.01, cv0=55.0, pv0=55.0, tau_planta=1.8, noise_pv=0.15,
)
L3 = LoopConfig(
    nombre="L3", kind="cascada", input_scans=40, dumax=1.0,
    orden=1, horizonte=5, sigma=20.0, delta=30.0, amp=20.0,
    ganancia_min=0.01, cv0=55.0, pv0=40.0, tau_planta=None, noise_pv=0.15,
    inner=F2, tau_tanque=70.0, noise_inner=1.5,
)
L4 = LoopConfig(
    nombre="L4", kind="cascada", input_scans=40, dumax=1.0,
    orden=1, horizonte=3, sigma=200.0, delta=30.0, amp=20.0,
    ganancia_min=0.01, cv0=55.0, pv0=40.0, tau_planta=None, noise_pv=0.15,
    inner=F1, tau_tanque=65.0, noise_inner=1.5,
)

LOOPS = {"F1": F1, "F2": F2, "L3": L3, "L4": L4}


# =============================================================================
# Planta: un paso de simulacion fisica (usada tanto en identificacion como
# en lazo cerrado). Devuelve el nuevo estado.
# =============================================================================
class Planta:
    """Modela: PID/HOREXT esclavo (si cascada) -> flujo (polo rapido) ->
    [estanque (polo lento), si cascada]. En directo, CV llega directo a flujo."""

    def __init__(self, cfg: LoopConfig, seed=0):
        self.cfg = cfg
        self.rng = np.random.default_rng(seed)
        self.flujo = cfg.inner.pv0 if cfg.kind == "cascada" else cfg.pv0
        self.nivel = cfg.pv0 if cfg.kind == "cascada" else None
        self.cv_interno = cfg.inner.cv0 if cfg.kind == "cascada" else cfg.cv0
        self.err_prev_interno = 0.0
        tau_flow = cfg.inner.tau_planta if cfg.kind == "cascada" else cfg.tau_planta
        self.a_flow = np.exp(-Tm / tau_flow)
        self.a_tank = np.exp(-Tm / cfg.tau_tanque) if cfg.kind == "cascada" else None
        self.cnt_inner = 0

    def step(self, cv_deseado, noise_on=True):
        """cv_deseado: para 'directo' es la salida del propio controlador
        (bomba). Para 'cascada' es el SP hacia el lazo interno (SP_F2/F1)."""
        cfg = self.cfg
        if cfg.kind == "directo":
            self.flujo = self.a_flow * self.flujo + (1 - self.a_flow) * cv_deseado
            pv = self.flujo + (self.rng.normal(0, cfg.noise_pv) if noise_on else 0.0)
            return pv, cv_deseado

        # ---- cascada: PID interno (siempre PID, igual que "usar HOREXT de
        # nivel con esclavo en PID", que es lo recomendado) ----
        inner = cfg.inner
        self.cnt_inner += 1
        Ts_inner = Tm * inner.input_scans
        if self.cnt_inner >= inner.input_scans:
            self.cnt_inner = 0
            fit_ruidoso = self.flujo + (self.rng.normal(0, cfg.noise_inner) if noise_on else 0.0)
            err = cv_deseado - fit_ruidoso
            rate_inner = inner.dumax * Ts_inner
            du = pid_du(1.6 if cfg.nombre == "L3" else 1.45, 0.75, Ts_inner,
                        err, self.err_prev_interno, rate_inner)
            self.err_prev_interno = err
            self.cv_interno = float(np.clip(self.cv_interno + du, 0, 95))

        self.flujo = self.a_flow * self.flujo + (1 - self.a_flow) * self.cv_interno
        self.nivel = self.a_tank * self.nivel + (1 - self.a_tank) * self.flujo
        pv = self.nivel + (self.rng.normal(0, cfg.noise_pv) if noise_on else 0.0)
        return pv, self.cv_interno


# =============================================================================
# MODO 1: identificacion en lazo abierto (PRBS + RLS)
# =============================================================================
def identificar(cfg: LoopConfig, T_end=None, noise_on=True, use_sensor_filter=False,
                 seed=1, quiet=False):
    if T_end is None:
        T_end = 900.0 if cfg.kind == "cascada" else 300.0

    prbs = PRBS(cfg.delta, cfg.amp, cfg.cv0)
    est = RLSHorext(cfg.orden, cfg.horizonte, cfg.sigma)
    planta = Planta(cfg, seed=seed)

    pv_prev = cfg.pv0
    cv_prev_lento = cfg.inner.cv0 if cfg.kind == "cascada" else cfg.cv0
    rate = cfg.dumax * Tm * cfg.input_scans

    cv_prev = cfg.cv0
    cv_man = cfg.cv0
    acumulador, n_acum = 0.0, 0
    cnt = 0
    n_steps = int(T_end / Tm)

    for k in range(n_steps):
        cv_man = prbs.step()

        du = float(np.clip(cv_man - cv_prev, -rate, rate))
        cv = cv_prev + du
        cv_prev = cv

        pv, cv_realizado = planta.step(cv, noise_on=noise_on)
        acumulador += pv
        n_acum += 1

        cnt += 1
        if cnt >= cfg.input_scans and k > 0:
            cnt = 0
            pv_usado = (acumulador / n_acum) if use_sensor_filter else pv
            acumulador, n_acum = 0.0, 0

            dy = pv_usado - pv_prev
            pv_prev = pv_usado
            du_realizado = cv_realizado - cv_prev_lento
            cv_prev_lento = cv_realizado

            est.update(dy, du_realizado)

    a0 = est.a[0]
    tau_impl = -(Tm * cfg.input_scans) / np.log(a0) if 0 < a0 < 1 else float("nan")

    if not quiet:
        etiqueta = "filtrado" if use_sensor_filter else "sin filtrar"
        print(f"[{cfg.nombre}] identificacion ({'con ruido' if noise_on else 'sin ruido'}, {etiqueta})")
        print(f"    a = {np.array2string(est.a, precision=3)}   "
              f"(a0 -> tau implicada = {tau_impl:6.1f} s)")
        print(f"    b = {np.array2string(est.b, precision=3)}")
        print(f"    lambda medio = {np.mean(est.lam_hist):.4f}   Traza P = {est.traza_p:.4f}")
        print()

    return dict(a=est.a, b=est.b, a0=a0, tau_impl=tau_impl,
                lam_medio=float(np.mean(est.lam_hist)), traza_p=est.traza_p)


# =============================================================================
# MODO 2: HOREXT en lazo cerrado con escalon de SP
# =============================================================================
def lazo_cerrado(cfg: LoopConfig, a, b, T_end=400.0, sp_step_time=60.0,
                   sp_step=-10.0, noise_on=True, adaptive=False, seed=1,
                   csv_path=None, quiet=False):
    """a, b: coeficientes ya identificados (arrays). Si adaptive=True, el
    estimador sigue corriendo en paralelo durante el lazo cerrado y el
    controlador usa el theta del CICLO LENTO ANTERIOR (replica el bug de
    orden JSR: LIC_XX corre antes que estimador_XX en MainRoutine)."""
    planta = Planta(cfg, seed=seed)
    rate = cfg.dumax * Tm * cfg.input_scans

    sp = cfg.pv0
    pv_prev = cfg.pv0
    cv_prev = cfg.inner.cv0 if cfg.kind == "cascada" else cfg.cv0
    cv = cv_prev

    est = RLSHorext(cfg.orden, cfg.horizonte, cfg.sigma) if adaptive else None
    cv_prev_lento = cv_prev
    a_usado, b_usado = np.array(a, dtype=float), np.array(b, dtype=float)

    t_log, pv_log, cv_log, sp_log, du_log, sat_log = [], [], [], [], [], []
    cnt = 0
    n_steps = int(T_end / Tm)

    for k in range(n_steps):
        t = k * Tm
        if abs(t - sp_step_time) < Tm / 2:
            sp = cfg.pv0 + sp_step

        pv, cv_realizado = planta.step(cv, noise_on=noise_on)

        cnt += 1
        if cnt >= cfg.input_scans and k > 0:
            cnt = 0
            dy = pv - pv_prev

            # --- controlador HOREXT usa el modelo "vigente" ---
            du, ganancia, a0_usado = horext_du(
                a_usado, b_usado, cfg.horizonte, dy, pv, sp,
                cfg.ganancia_min, rate)
            cv = float(np.clip(cv_prev + du, 0, 95))
            sat = abs(du) >= rate * 0.999

            # --- estimador adaptativo (opcional), corre DESPUES del
            # controlador dentro del mismo scan -> replica el bug de orden ---
            if adaptive:
                du_realizado = cv_realizado - cv_prev_lento
                cv_prev_lento = cv_realizado
                est.update(dy, du_realizado)
                a_usado, b_usado = est.a, est.b   # el PROXIMO ciclo usara esto

            pv_prev = pv
            cv_prev = cv

            t_log.append(t); pv_log.append(pv); cv_log.append(cv)
            sp_log.append(sp); du_log.append(du); sat_log.append(sat)

    pv_arr = np.array(pv_log); sp_arr = np.array(sp_log)
    t_arr = np.array(t_log)
    idx_step = np.searchsorted(t_arr, sp_step_time)
    post = pv_arr[idx_step:]
    sp_post = sp_arr[idx_step:]
    err_post = sp_post - post

    if sp_step < 0:
        overshoot = max(0.0, sp_post[0] - post.min())   # bajando: undershoot = pasarse hacia abajo
    else:
        overshoot = max(0.0, post.max() - sp_post[0])

    banda = 0.02 * abs(sp_step)
    settled_mask = np.abs(err_post) < max(banda, 1e-6)
    t_settle = float("nan")
    for i in range(len(settled_mask)):
        if settled_mask[i] and np.all(settled_mask[i:]):
            t_settle = t_arr[idx_step + i] - sp_step_time
            break

    # cruces "reales": ignora el dither de ruido del sensor cerca de cero,
    # solo cuenta vueltas que superen 3-sigma del ruido propio del PV
    tail = err_post[len(err_post) // 2:]
    epsilon = max(banda, 3.0 * cfg.noise_pv)
    signo = np.where(np.abs(tail) > epsilon, np.sign(tail), 0)
    signo_no_cero = signo[signo != 0]
    cruces = int(np.sum(np.diff(signo_no_cero) != 0)) if len(signo_no_cero) > 2 else 0
    sat_frac = float(np.mean(sat_log[idx_step:])) if len(sat_log) > idx_step else 0.0

    resultado = dict(overshoot=float(overshoot), t_settle=t_settle,
                      cruces_oscilacion=cruces, sat_frac=sat_frac)

    if not quiet:
        print(f"[{cfg.nombre}] lazo cerrado HOREXT (escalon {sp_step:+.1f}%, "
              f"{'adaptativo' if adaptive else 'modelo fijo'}, "
              f"{'con ruido' if noise_on else 'sin ruido'})")
        tipo = "undershoot" if sp_step < 0 else "overshoot"
        print(f"    {tipo} = {overshoot:.2f} %   t_asentamiento = {t_settle:.1f} s")
        print(f"    cruces de error tras el transiente = {cruces}  "
              f"({'posible oscilacion sostenida' if cruces >= 3 else 'ok'})")
        print(f"    saturado en rate-limit {sat_frac*100:.0f}% del tiempo post-escalon")
        print()

    if csv_path:
        import csv
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["t", "SP", "PV", "CV", "du", "saturado"])
            w.writerows(zip(t_log, sp_log, pv_log, cv_log, du_log, sat_log))

    return resultado


# =============================================================================
# Bateria: corre identificacion + lazo cerrado para los 4 lazos
# =============================================================================
def bateria():
    for nombre, cfg in LOOPS.items():
        print("=" * 74)
        print(f"LAZO {nombre}")
        print("=" * 74)
        r_sim = identificar(cfg, noise_on=False)
        r_real_sf = identificar(cfg, noise_on=True, use_sensor_filter=False)
        r_real_f = identificar(cfg, noise_on=True, use_sensor_filter=True)

        # usa el resultado con filtro (mas confiable) para el lazo cerrado
        paso = -10.0 if cfg.kind == "cascada" else 10.0
        lazo_cerrado(cfg, r_real_f["a"], r_real_f["b"], sp_step=paso, noise_on=True)


# =============================================================================
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="modo", required=True)

    p1 = sub.add_parser("ident", help="identificacion en lazo abierto")
    p1.add_argument("--loop", choices=LOOPS.keys(), required=True)
    p1.add_argument("--no-noise", action="store_true")
    p1.add_argument("--filtro", action="store_true", help="promediar sensor en el ciclo lento")
    p1.add_argument("--seed", type=int, default=1)

    p2 = sub.add_parser("step", help="HOREXT en lazo cerrado con escalon de SP")
    p2.add_argument("--loop", choices=LOOPS.keys(), required=True)
    p2.add_argument("--sp-step", type=float, default=-10.0)
    p2.add_argument("--no-noise", action="store_true")
    p2.add_argument("--adaptive", action="store_true", help="estimador vivo durante el lazo cerrado (replica bug JSR)")
    p2.add_argument("--csv", type=str, default=None)
    p2.add_argument("--seed", type=int, default=1)

    sub.add_parser("bateria", help="corre todo, todos los lazos")

    args = ap.parse_args()

    if args.modo == "ident":
        identificar(LOOPS[args.loop], noise_on=not args.no_noise,
                     use_sensor_filter=args.filtro, seed=args.seed)
    elif args.modo == "step":
        cfg = LOOPS[args.loop]
        r = identificar(cfg, noise_on=not args.no_noise, use_sensor_filter=True,
                          seed=args.seed, quiet=True)
        lazo_cerrado(cfg, r["a"], r["b"], sp_step=args.sp_step,
                       noise_on=not args.no_noise, adaptive=args.adaptive,
                       csv_path=args.csv, seed=args.seed)
    elif args.modo == "bateria":
        bateria()
