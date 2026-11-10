# 🧪 Real Gas Thermodynamics Simulator (Peng-Robinson EOS)

An interactive thermodynamic property calculator and isotherm visualizer built in Python using **Streamlit**, **NumPy**, and **Plotly**. This software models the non-ideal behavior of pure real gases across wide temperature and pressure ranges using the **Peng-Robinson (1976) Cubic Equation of State**.

🔗 **Live Deployed Web Application:** [https://real-gas-thermo-raghav.streamlit.app/](https://real-gas-thermo-raghav.streamlit.app/)

---

## 📌 Project Overview & Modules

This software was developed across three core thermodynamic modules:

1. **Module 1 — Base Software (PVT Behavior):** Solves the cubic polynomial for the Compressibility Factor ($Z$) and computes the real gas Molar Volume ($V$) vs. Ideal Gas Molar Volume ($V^{ig}$).
2. **Module 2 — Homework Extension (Thermodynamic Departure Properties):** Evaluates residual Enthalpy ($H - H^{ig}$) and residual Entropy ($S - S^{ig}$) using analytical temperature derivatives of the attractive parameter $a(T)$.
3. **Module 3 — Computer Project (Pure Component Fugacity):** Determines the dimensionless fugacity coefficient ($\phi$) and pure component fugacity ($f$) to quantify chemical potential deviation from ideal gas behavior, accompanied by interactive $Z$ vs. $P$ and $\phi$ vs. $P$ isotherm plots.

---

## 📐 Mathematical Formulation

### 1. Peng-Robinson Equation of State
The explicit pressure form of the Peng-Robinson EOS is given by:

$$P = \frac{RT}{V - b} - \frac{a(T)}{V(V + b) + b(V - b)}$$

Where the pure-component critical parameters and temperature-dependent $\alpha$-function are defined as:

$$a_c = 0.45723553 \frac{R^2 T_c^2}{P_c}, \quad b = 0.07779607 \frac{R T_c}{P_c}$$

$$\alpha(T_r, \omega) = \left[1 + \kappa \left(1 - \sqrt{T_r}\right)\right]^2, \quad T_r = \frac{T}{T_c}$$

$$\kappa = 0.37464 + 1.54226\omega - 0.26992\omega^2$$

$$a(T) = a_c \cdot \alpha(T_r, \omega)$$

The analytical temperature derivative $\frac{da}{dT}$ required for Enthalpy and Entropy calculations is:

$$\frac{da}{dT} = -a_c \kappa \sqrt{\frac{\alpha}{T \cdot T_c}}$$

---

### 2. Polynomial Form for Compressibility Factor ($Z$)
Defining dimensionless parameters $A$ and $B$:

$$A = \frac{a(T)P}{R^2 T^2}, \quad B = \frac{bP}{RT}$$

The EOS is expressed as a cubic polynomial in $Z = \frac{PV}{RT}$:

$$Z^3 - (1 - B)Z^2 + (A - 2B - 3B^2)Z - (AB - B^2 - B^3) = 0$$

For gas/supercritical vapor phase modeling, the largest real root satisfying $Z > B$ is extracted using `numpy.roots`, and the real molar volume is obtained via $V = \frac{ZRT}{P}$.

---

### 3. Enthalpy and Entropy Departure Functions
Departure (residual) properties measure the difference between the real fluid state $(T, P)$ and the ideal gas state at the same $(T, P)$:

$$\frac{H - H^{ig}}{RT} = (Z - 1) + \frac{T\left(\frac{da}{dT}\right) - a}{2\sqrt{2}bRT} \ln\left[\frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B}\right]$$

$$\frac{S - S^{ig}}{R} = \ln(Z - B) + \frac{\frac{da}{dT}}{2\sqrt{2}bR} \ln\left[\frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B}\right]$$

---

### 4. Pure Component Fugacity ($f$)
The fugacity coefficient ($\phi = \frac{f}{P}$) is computed from the isothermal residual Gibbs free energy relation:

$$\ln \phi = (Z - 1) - \ln(Z - B) - \frac{A}{2\sqrt{2}B} \ln\left[\frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B}\right]$$

$$f = \phi \cdot P$$

---

## 🗂️ Preset Pure Component Database

The application includes built-in critical constants for six industrial fluids alongside a **Custom Gas** mode:

| Fluid | Formula | $T_c$ (K) | $P_c$ (bar) | Acentric Factor ($\omega$) |
| :--- | :--- | :---: | :---: | :---: |
| **Carbon Dioxide** | $\text{CO}_2$ | 304.12 | 73.74 | 0.224 |
| **Methane** | $\text{CH}_4$ | 190.56 | 45.99 | 0.011 |
| **Nitrogen** | $\text{N}_2$ | 126.20 | 33.98 | 0.037 |
| **Ethane** | $\text{C}_2\text{H}_6$ | 305.32 | 48.72 | 0.099 |
| **Propane** | $\text{C}_3\text{H}_8$ | 369.83 | 42.48 | 0.152 |
| **Water / Steam** | $\text{H}_2\text{O}$ | 647.14 | 220.64 | 0.344 |

---
