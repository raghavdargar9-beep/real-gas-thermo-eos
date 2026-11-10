# 🧪 Unified Cubic EOS Simulator: PVT, H, S & Fugacity of Real Gases

**Author:** Raghav Dargar (Entry No: 2025CH11790)  
🔗 **Live Deployed Web Application:** [https://real-gas-thermo-raghav.streamlit.app/](https://real-gas-thermo-raghav.streamlit.app/)

An interactive chemical engineering thermodynamics web application built in Python using **Streamlit**, **NumPy**, and **Plotly**. This software predicts the PVT behavior, residual and absolute Enthalpy ($H$) and Entropy ($S$), pure-component fugacity ($f$), and vapor-liquid phase stability of real gases using four unified Cubic Equations of State:
* **van der Waals (vdW)**
* **Redlich-Kwong (RK)**
* **Soave-Redlich-Kwong (SRK)**
* **Peng-Robinson (PR)**

---

## 📌 Project Modules & Interactive Visualizations

1. **Part 1 — Base Software (PVT Behavior):** Solves the generalized cubic compressibility polynomial analytically via the **Cardano-Vieta trigonometric method** to extract single-phase or two-phase ($Z_{liquid}, Z_{vapor}$) roots and molar volumes ($V_m$).
2. **Part 2 — Homework Extension (Thermodynamic Properties $H$ & $S$):** Computes residual enthalpy ($H^R$) and residual entropy ($S^R$) via density integrals, and integrates ideal-gas heat capacity polynomials ($C_p^{ig}/R$) from $T_{ref} = 298.15\text{ K}, P_{ref} = 1\text{ bar}$ to determine absolute $H$ and $S$.
3. **Part 3 — Computer Project (Pure Component Fugacity & Phase Stability):** Evaluates the fugacity coefficient ($\phi$), pure component fugacity ($f = \phi P$), and normalized residual Gibbs energy ($G^R/RT = \ln\phi$) to automatically identify the thermodynamically stable phase.
4. **Interactive Plots Included in the Web App:**
   * **Log-Scale $P$-$V_m$ Isotherm Plot:** Tracks real-fluid compression and highlights the operating state.
   * **Pure Gas Fugacity ($f$) vs. System Pressure ($P$) Curve:** Visualizes chemical potential departure from the ideal gas reference line ($f = P$).
   * **Multi-EOS Comparison Plots ($Z$ vs. $P$ and $H^R$ vs. $P$):** Overlays `vdW`, `RK`, `SRK`, and `PR` isotherms on a single canvas to illustrate model divergence at high pressure.

---

## 📐 Unified Mathematical Formulation

### 1. Generalized Cubic Equation of State
All four models are expressed in the unified form (Smith, Van Ness, Abbott — *Introduction to Chemical Engineering Thermodynamics*, 8th Ed.):

$$P = \frac{RT}{V - b} - \frac{a(T)}{(V + \epsilon b)(V + \sigma b)}$$

$$a(T) = \Psi \frac{\alpha(T_r) R^2 T_c^2}{P_c}, \quad b = \Omega \frac{R T_c}{P_c}$$

| EOS Model | $\alpha(T_r)$ | $\sigma$ | $\epsilon$ | $\Omega$ | $\Psi$ | $Z_c$ (Predicted) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **vdW** | $1$ | $0$ | $0$ | $1/8$ | $27/64$ | $0.375$ |
| **RK** | $T_r^{-1/2}$ | $1$ | $0$ | $0.08664$ | $0.42748$ | $0.333$ |
| **SRK** | $[1 + \kappa_{SRK}(1 - \sqrt{T_r})]^2$ | $1$ | $0$ | $0.08664$ | $0.42748$ | $0.333$ |
| **PR** | $[1 + \kappa_{PR}(1 - \sqrt{T_r})]^2$ | $1 + \sqrt{2}$ | $1 - \sqrt{2}$ | $0.07780$ | $0.45724$ | $0.307$ |

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

## ✅ Textbook Benchmark Validation

### 1. SVNA Example 3.9 (n-Butane via RK EOS at $350\text{ K}, 9.4573\text{ bar}$)
* **Vapor Molar Volume ($V_v$):** Program predicts $2555.3\text{ cm}^3\text{/mol}$ vs. SVNA textbook value of $2555\text{ cm}^3\text{/mol}$.
* **Liquid Molar Volume ($V_l$):** Program predicts $133.3\text{ cm}^3\text{/mol}$ vs. SVNA textbook value of $133.3\text{ cm}^3\text{/mol}$.

### 2. SVNA Table 13.7 (n-Butane at $500\text{ K}, 50\text{ bar}$)

| EOS Model | $Z$ (SVNA Book) | $Z$ (Program) | $H^R$ Book (J/mol) | $H^R$ Program (J/mol) | $S^R$ Book (J/mol·K) | $S^R$ Program (J/mol·K) | $\phi$ (Program) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **vdW** | 0.6608 | 0.6610 | -3937 | -3935.0 | -5.66 | -5.654 | 0.6788 |
| **RK** | 0.6850 | 0.6852 | -4505 | -4503.0 | -6.58 | -6.576 | 0.6773 |
| **SRK** | 0.7222 | 0.7224 | -4824 | -4821.5 | -7.24 | -7.239 | 0.6989 |
| **PR** | 0.6907 | 0.6909 | -4988 | -4985.2 | -7.43 | -7.426 | 0.6708 |

---

## ⚠️ Discussion of Limitations & Concrete Examples

While SRK and PR are industry standards for non-polar hydrocarbons (predicting vapor pressures within $0.3\%$ for Propane at $300\text{ K}$ and $\text{CO}_2$ at $280\text{ K}$, and $H^R$ within $\sim 100\text{ J/mol}$), two-parameter cubic Equations of State exhibit four fundamental physical limitations:

1. **Overprediction of Critical Compressibility ($Z_c$):**
   * **Theory:** Two-parameter cubic EOS models fix $a_c$ and $b$ solely from $T_c$ and $P_c$, forcing a universal constant critical compressibility factor $Z_c = \frac{P_c V_c}{R T_c}$ (`vdW` = $0.375$, `RK`/`SRK` = $0.333$, `PR` = $0.307$).
   * **Concrete Example:** Real fluids have much lower experimental $Z_c$ values. For **n-Butane** ($Z_c^{exp} = 0.274$) and **$\text{CO}_2$** ($Z_c^{exp} = 0.274$), PR overpredicts the critical molar volume by **$+12.0\%$** and vdW overpredicts by **$+36.9\%$**. For strongly polar **Water** ($Z_c^{exp} = 0.229$), PR overpredicts critical volume by **$+34.1\%$** and vdW by **$+63.8\%$**.

2. **Inadequacy for Polar & Hydrogen-Bonding Fluids (Water, Alcohols, Ammonia):**
   * **Theory:** Standard cubic EOS models assume spherical or weakly non-spherical van der Waals dispersion forces parameterized solely by reduced temperature ($T_r$) and Pitzer's acentric factor ($\omega$). They lack association terms for directional hydrogen bonding or permanent dipole-dipole interactions.
   * **Concrete Example:** For **Water ($\text{H}_2\text{O}$) at $373.15\text{ K}$ ($100^\circ\text{C}$)**, the experimental saturation vapor pressure is **$1.013\text{ bar}$**, whereas standard Peng-Robinson predicts **$0.960\text{ bar}$ ($5.3\%$ error)**. Furthermore, at $373.15\text{ K}$ and $1.013\text{ bar}$, PR predicts a saturated liquid molar volume of **$V_l = 22.4\text{ cm}^3\text{/mol}$** compared to the experimental steam-table value of **$18.8\text{ cm}^3\text{/mol}$ ($+19.1\%$ error)**.

3. **Failure of Redlich-Kwong (RK) for Vapor Pressures & Liquid Volumes:**
   * **Theory:** Because RK uses $\alpha(T_r) = T_r^{-1/2}$ with zero dependence on the acentric factor ($\omega$), it cannot account for molecular chain length or non-sphericity in the attractive term.
   * **Concrete Example:** While RK achieves low error for simple gas-phase compressibility ($Z = 0.6852$ for n-Butane at $500\text{ K}, 50\text{ bar}$), its subcritical vapor pressure predictions for heavier hydrocarbons ($C_4+$) deviate by **$>15\%$**, whereas SRK and PR reduce this error to **$<0.5\%$** by introducing $\kappa(\omega)$.

4. **Temperature Extrapolation Limit of Ideal-Gas Heat Capacity ($C_p^{ig}$) Polynomials:**
   * **Theory:** The polynomial correlation $\frac{C_p^{ig}}{R} = A + BT + CT^2 + DT^{-2}$ from SVNA Table C.1 is empirically regressed only from **$298.15\text{ K}$ to $T_{max}$ ($1500\text{ K}$)**.
   * **Concrete Example:** For fluids with a non-zero $D$ constant such as **$\text{CO}_2$** ($D = -1.157 \times 10^5\text{ K}^2$), extrapolating the absolute enthalpy/entropy integral into the cryogenic region (e.g., $T = 150\text{ K}$) causes the $D T^{-2}$ term to quadruple in magnitude, producing unphysical heat capacities.
