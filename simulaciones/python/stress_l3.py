"""
Stress test of the L3 identification path as implemented in 4tanks.L5K.

Replicates, at the 0.1s scan level:
  - PRBS_L3 (7-bit LFSR, taps 6,7; seed 0011100; delta=30s; amp=20)
  - LIC_03 manual branch: du_L3 = CV_L3_man - CV_L3_prev, clamp +-rate (dumax_L3=1%/s * 4s)
  - selector_2 with CV_L2=0, Sel2=1: CV_L3_prev := FIT_02 EVERY SCAN  <-- the stomp
  - FIC_02 PID (Kc=1.6, Ki=0.75, Ts=1s) with dumax_F2=1%/s rate limiter
  - Plant: flow loop ~ first order + tank integrator/first-order (tau configurable)
  - estimador_L3: RLS on increments, orden=1, horizonte=5 (n_beta=5), Ts=4s,
    du input = du_F2_lento (realized pump increments), Fortescue lambda, Sigma=20
Scenarios:
  A) sim-like: no sensor noise, stomp ON   (what the user validated)
  B) real-like: flow+level noise, stomp ON (what the plant does)
  C) real-like: flow+level noise, stomp OFF (selector tracking fixed)
Also reports the effective PRBS amplitude actually reaching SP_F2.
"""
import numpy as np

rng = np.random.default_rng(1)

def run(scenario, T_end=600.0, tau_tank=70.0, tau_flow=8.0,
        noise_flow=0.0, noise_level=0.0, stomp=True, seed=1):
    rng = np.random.default_rng(seed)
    dt = 0.1
    n_steps = int(T_end/dt)

    # PRBS state
    l = [0,0,1,1,1,0,0]
    act = 0
    deltaL3, amp = 30.0, 20.0
    Tm = 0.1
    CV_L3_man = 55.0
    pmax, pmin = CV_L3_man+amp/2, CV_L3_man-amp/2
    PL3_out = CV_L3_man

    # LIC_03 manual
    CV_L3_prev = 55.0
    CV_L3 = 55.0
    dumax_L3 = 1.0
    input_scans_L3 = 40
    rate_L3 = dumax_L3 * Tm*input_scans_L3   # 4 % per slow cycle

    # FIC_02 PID
    Kc, Ki = 1.6, 0.75
    Ts_F2 = 1.0
    dumax_F2 = 1.0
    rate_F2 = dumax_F2 * Ts_F2
    CV_F2 = 55.0
    e_prev = 0.0
    SP_F2 = 55.0

    # plant: pump CV -> flow (first order tau_flow), flow -> level
    # level (linearized): tau_tank, steady-state gain 1 (percent to percent)
    flow = 55.0
    level = 40.0
    a_f = np.exp(-dt/tau_flow)
    a_t = np.exp(-dt/tau_tank)

    # estimator
    orden, horizonte = 1, 5
    n_beta = orden + horizonte - 1
    npar = 2*orden + horizonte - 1
    theta = np.zeros(npar)
    P = np.eye(npar)*10000.0
    dy_hist = np.zeros(10)
    du_hist = np.zeros(10)
    LIT_prev = level
    CVF2_prev_lento = CV_F2
    Sigma = 20.0
    lam_hist = []
    a0_hist = []

    cnt_L3 = 0
    cnt_F2 = 0
    sp_f2_log = []

    for k in range(n_steps):
        # ---- PRBS every scan
        act -= 1
        if act <= 0:
            act = int(deltaL3/Tm)
            s = l[5] ^ l[6]
            l = [s] + l[:6]
            PL3_out = pmax if s else pmin
        CV_L3_man = PL3_out

        # ---- selector_2 stomp (runs every scan, CV_L2=0 always loses... wins)
        FIT_02 = flow + (rng.normal(0, noise_flow) if noise_flow else 0.0)
        if stomp:
            CV_L3_prev = FIT_02   # CV_L2=0 <= CV_L3 -> this branch every scan

        # ---- LIC_03 manual, slow cycle (every 40 scans)
        cnt_L3 += 1
        if cnt_L3 >= input_scans_L3:
            cnt_L3 = 0
            du = CV_L3_man - CV_L3_prev
            du = np.clip(du, -rate_L3, rate_L3)
            CV_L3 = CV_L3_prev + du
            CV_L3_prev = CV_L3
            sp_f2_log.append(CV_L3)

        SP_F2 = CV_L3

        # ---- FIC_02 PID every 10 scans (Ts=1s)
        cnt_F2 += 1
        if cnt_F2 >= 10:
            cnt_F2 = 0
            e = SP_F2 - FIT_02
            du2 = Kc*(1+Ki*Ts_F2/2)*e - Kc*(1-Ki*Ts_F2/2)*e_prev
            e_prev = e
            du2 = np.clip(du2, -rate_F2, rate_F2)
            CV_F2 = np.clip(CV_F2 + du2, 0, 95)

        # ---- plant
        flow = a_f*flow + (1-a_f)*CV_F2
        level = a_t*level + (1-a_t)*flow

        # ---- estimador_L3 slow cycle (fires same scans as LIC_03; runs after)
        if cnt_L3 == 0 and k > 0:
            LIT_03 = level + (rng.normal(0, noise_level) if noise_level else 0.0)
            dy = LIT_03 - LIT_prev
            LIT_prev = LIT_03
            du_l = CV_F2 - CVF2_prev_lento
            CVF2_prev_lento = CV_F2

            phi = np.zeros(npar)
            phi[0] = -dy_hist[0]
            phi[1:1+n_beta] = du_hist[:n_beta]

            est = phi @ theta
            e_pred = dy - est
            Pphi = P @ phi
            lam_dem = 1.0 + phi @ Pphi
            lam = 1.0 - e_pred**2/(lam_dem*Sigma)
            lam = min(max(lam, 0.5), 1.0)
            den = lam + phi @ Pphi
            if den > 1e-6:
                K = Pphi/den
                theta = theta + K*e_pred
                P = (P - np.outer(K, Pphi))/lam
            dy_hist[1:] = dy_hist[:-1]; dy_hist[0] = dy
            du_hist[1:] = du_hist[:-1]; du_hist[0] = du_l
            lam_hist.append(lam)
            a0_hist.append(-theta[0])

    a0 = -theta[0]
    bs = theta[1:1+n_beta]
    tau_implied = -4.0/np.log(a0) if 0 < a0 < 1 else float('nan')
    sp = np.array(sp_f2_log)
    print(f"{scenario}")
    print(f"  a0 = {a0:+.3f}  (tau implicada = {tau_implied:5.1f} s ; tau real tanque = {tau_tank} s)")
    print(f"  b  = {np.array2string(bs, precision=3)}")
    print(f"  lambda medio = {np.mean(lam_hist):.4f}   TrazaP = {np.trace(P):.3f}")
    print(f"  SP_F2 efectivo: min {sp.min():.1f}  max {sp.max():.1f}  (PRBS pedia {pmin:.0f} a {pmax:.0f})")
    print()

print("="*74)
print("A) SIMULACION (sin ruido, stomp selector ACTIVO) - lo que valida el user")
run("A", noise_flow=0.0, noise_level=0.0, stomp=True)

print("B) PLANTA REAL (ruido flujo sigma=1.5%, nivel sigma=0.15%, stomp ACTIVO)")
run("B", noise_flow=1.5, noise_level=0.15, stomp=True)

print("C) PLANTA REAL, MISMO RUIDO, pero SIN el stomp del selector (corregido)")
run("C", noise_flow=1.5, noise_level=0.15, stomp=False)

print("D) PLANTA REAL sin stomp + flujo mas lento (tau_flow=12s, rate-limit domina)")
run("D", noise_flow=1.5, noise_level=0.15, stomp=False, tau_flow=12.0)
