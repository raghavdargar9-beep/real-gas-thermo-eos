# 🧪 Unified Cubic EOS Simulator: PVT, H, S & Fugacity of Real Gases

**Author:** Raghav Dargar (Entry No: 2025CH11790)  
🔗 **Live Deployed Web Application:** [https://real-gas-thermo-raghav.streamlit.app/](https://real-gas-thermo-raghav.streamlit.app/)

An interactive chemical engineering thermodynamics web application built in Python using **Streamlit**, **NumPy**, and **Plotly**. This software predicts the PVT behavior, residual and absolute Enthalpy ($H$) and Entropy ($S$), pure-component fugacity ($f$), and vapor-liquid phase stability of real gases using four unified Cubic Equations of State:
* **van der Waals (vdW)**
* **Redlich-Kwong (RK)**
* **Soave-Redlich-Kwong (SRK)**
* **Peng-Robinson (PR)**

---

## 📌 Project Modules

1. **Part 1 — Base Software (PVT Behavior):** Solves the generalized cubic compressibility polynomial analytically via the **Cardano-Vieta trigonometric method** to extract single-phase or two-phase ($Z_{liquid}, Z_{vapor}$) roots and molar volumes ($V_m$).
2. **Part 2 — Homework Extension (Thermodynamic Properties $H$ & $S$):** Computes residual enthalpy ($H^R$) and residual entropy ($S^R$) via density integrals, and integrates ideal-gas heat capacity polynomials ($C_p^{ig}/R$) from $T_{ref} = 298.15\text{ K}, P_{ref} = 1\text{ bar}$ to determine absolute $H$ and $S$.
3. **Part 3 — Computer Project (Pure Component Fugacity & Phase Stability):** Evaluates the fugacity coefficient ($\phi$), pure component fugacity ($f = \phi P$), and normalized residual Gibbs energy ($G^R/RT = \ln\phi$) to automatically identify the thermodynamically stable phase.

---

## 📐 Unified Mathematical Formulation

### 1. Generalized Cubic Equation of State
All four models are expressed in the unified form (Smith, Van Ness, Abbott — *Introduction to Chemical Engineering Thermodynamics*, 8th Ed.):

$$P = \frac{RT}{V - b} - \frac{a(T)}{(V + \epsilon b)(V + \sigma b)}$$

$$a(T) = \Psi \frac{\alpha(T_r) R^2 T_c^2}{P_c}, \quad b = \Omega \frac{R T_c}{P_c}$$

| EOS Model | $\alpha(T_r)$ | $\sigma$ | $\epsilon$ | $\Omega$ | $\Psi$ |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **vdW** | $1$ | $0$ | $0$ | $1/8$ | $27/64$ |
| **RK** | $T_r^{-1/2}$ | $1$ | $0$ | $0.08664$ | $0.42748$ |
| **SRK** | $[1 + \kappa_{SRK}(1 - \sqrt{T_r})]^2$ | $1$ | $0$ | $0.08664$ | $0.42748$ |
| **PR** | $[1 + \kappa_{PR}(1 - \sqrt{T_r})]^2$ | $1 + \sqrt{2}$ | $1 - \sqrt{2}$ | $0.07780$ | $0.45724$ |

---

### 2. Compressibility Polynomial & Analytical Cardano-Vieta Solver
Defining $A = \frac{aP}{(RT)^2}$, $B = \frac{bP}{RT}$, $u = \sigma + \epsilon$, and $w = \sigma\epsilon$:

$$Z^3 + (-1 - B + uB)Z^2 + (A + wB^2 - uB - uB^2)Z + (-AB - wB^2 - wB^3) = 0$$

Physical roots satisfy $Z > B$. When three real roots exist in the subcritical region, the program evaluates $\frac{G^R}{RT} = \ln\phi$ for both the liquid ($Z_{min}$) and vapor ($Z_{max}$) roots and selects the root yielding the **minimum Gibbs energy**.

---

### 3. Residual Properties, Absolute $H$ & $S$, and Fugacity
Using the generalized integral $I = \frac{1}{\sigma - \epsilon}\ln\left(\frac{Z + \sigma B}{Z + \epsilon B}\right)$ (and $I = \frac{B}{Z}$ for vdW) with $q = \frac{A}{B}$:

$$\frac{H^R}{RT} = (Z - 1) + \left[\frac{d\ln\alpha}{d\ln T_r} - 1\right] q I$$

$$\frac{S^R}{R} = \ln(Z - B) + \left[\frac{d\ln\alpha}{d\ln T_r}\right] q I$$

$$\ln\phi = \frac{G^R}{RT} = (Z - 1) - \ln(Z - B) - q I \implies f = \phi P$$

Absolute enthalpy and entropy relative to the ideal gas reference state ($T_{ref} = 298.15\text{ K}, P_{ref} = 1.0\text{ bar}$) are computed using Table C.1 polynomial constants ($\frac{C_p^{ig}}{R} = A + BT + CT^2 + DT^{-2}$):

$$H(T, P) = R\int_{298.15}^{T} \frac{C_p^{ig}}{R} dT + H^R(T, P)$$

$$S(T, P) = R\int_{298.15}^{T} \frac{C_p^{ig}}{R} \frac{dT}{T} - R\ln\left(\frac{P}{1.0}\right) + S^R(T, P)$$

---

## ✅ Textbook Benchmark Validation (SVNA Table 13.7: n-Butane at 500 K, 50 bar)

| EOS Model | $Z$ (SVNA Book) | $Z$ (Program) | $H^R$ Book (J/mol) | $H^R$ Program (J/mol) | $\phi$ (Program) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **vdW** | 0.6608 | 0.6610 | -3937 | -3935.0 | 0.6788 |
| **RK** | 0.6850 | 0.6852 | -4505 | -4503.0 | 0.6773 |
| **SRK** | 0.7222 | 0.7224 | -4824 | -4821.5 | 0.6989 |
| **PR** | 0.6907 | 0.6909 | -4988 | -4985.2 | 0.6708 |
