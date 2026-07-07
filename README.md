# Control System for a 4-Tank Pilot Plant

![Status](https://img.shields.io/badge/Status-Completed-brightgreen)
![University](https://img.shields.io/badge/UdeC-Electronic_Civil_Engineering-orange)
![Course](https://img.shields.io/badge/Course-547507_Computer_Process_Control-lightgrey)

🇪🇸 [Versión en Español](README.es.md)

## Overview
This repository documents the design, simulation, and implementation of an automatic control system for a 4-tank pilot plant. The project spans mathematical modeling and discretization of continuous systems, through the synthesis of advanced control algorithms, all the way to direct deployment on industrial hardware.

▶️ **Demo Videos:** 
[Watch the PID Control loop + Feedforward applied on the real plant](https://youtu.be/Ni7OxD4VFhw)
[Watch the PRBS identification and closed-loop HOREXT control running on the real plant](https://youtu.be/QE6KxWTnEUY)

📄 **Full Technical Report:** [Technical Report (PDF)](docs/documentacion/main.pdf) — mathematical modeling, controller synthesis, and a full real-time troubleshooting log.

### MIMO Operating-Point Simulator
> [Interactive simulator — click here](https://aguscsc.github.io/Sistema-de-Control-para-Planta-Piloto-de-4-Estanques/simulaciones/index.html)

## HMI
![HMI](pics/planta.gif)

*Human-Machine Interface designed under the ISA-101 standard.*

## Control Architecture & P&ID
![P&ID](pics/instrumentacion/planta_4_estanques_P&ID.png)
*Instrumentation diagram based on ISA-5.1, detailing the level (LIC) and flow (FIC) control loops.*

## Technical Scope
* **Real-Time Control:** Design and tuning of discrete single-loop controllers for the tank cascade dynamics.
* **Adaptive Predictive Control:** Extended-Horizon Self-Tuning Control (HOREXT), derived independently from Ydstie, Kershenbaum & Sargent (1985) — online RLS identification with Fortescue's variable forgetting factor, extrapolated to a multi-step predictor, tuned and validated in closed loop on the physical plant.
* **Industrial Communications:** Data transmission protocols for robust synchronization between instrumentation, PLC, and computing equipment.
* **HMI Design:** Graphical interfaces for operation and monitoring, structured under the high-performance **ISA-101** standard.

## Tools & Technologies
* **Control Hardware:** Allen-Bradley ControlLogix 1756-L81E.
* **PLC Programming:** Studio 5000 (Rockwell Automation), Structured Text (IEC 61131-3).
* **HMI Development:** FactoryTalk (Rockwell Automation).
* **Communications:** EtherNet/IP and an OPC UA server (FactoryTalk Gateway) for integration with MATLAB/Simulink.
* **Simulation & Analysis:** MATLAB / Simulink, Python.

---
## Engineers
- **[Agustín Torres](https://github.com/aguscsc)**
- **[Ignacio Cerda](https://github.com/LovesCharlie)**
- **[Leví Sojos](https://github.com/gadivalr)**

---
