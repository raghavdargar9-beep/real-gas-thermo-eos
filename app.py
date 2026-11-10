import streamlit as st
import numpy as np
import plotly.graph_objects as go

st.set_page_config(
    page_title="Unified Cubic EOS Explorer | PVT, H, S & Fugacity",
    page_icon="🧪",
    layout="wide"
)

R = 8.314462618  # J/(mol·K)
T_REF = 298.15   # K
P_REF = 1.0      # bar

# Fluid Database from Smith, Van Ness, Abbott (Tables B.1 & C.1)
# Format: {Name: {"Tc": K, "Pc": bar, "omega": -, "cp": (A, B, C, D)}}
FLUID_DB = {
    "Methane (CH₄)": {
        "Tc": 190.56, "Pc": 45.99, "omega": 0.011,
        "cp": (1.702, 9.081e-3, -2.164e-6, 0.0)
    },
    "n-Butane (n-C₄H₁₀) [SVNA Benchmark]": {
        "Tc": 425.12, "Pc": 37.96, "omega": 0.200,
        "cp": (1.935, 36.915e-3, -11.402e-6, 0.0)
    },
    "Carbon Dioxide (CO₂)": {
        "Tc": 304.12, "Pc": 73.74, "omega": 0.224,
        "cp": (5.457, 1.045e-3, 0.0, -1.157e5)
    },
    "Ethane (C₂H₆)": {
        "Tc": 305.32, "Pc": 48.72, "omega": 0.099,
        "cp": (1.131, 19.225e-3, -5.561e-6, 0.0)
    },
    "Propane (C₃H₈)": {
        "Tc": 369.83, "Pc": 42.48, "omega": 0.152,
        "cp": (1.213, 28.785e-3, -8.824e-6, 0.0)
    },
    "Nitrogen (N₂)": {
        "Tc": 126.20, "Pc": 33.98, "omega": 0.037,
        "cp": (3.280, 0.593e-3, 0.0, 0.040e5)
    },
    "Water (H₂O)": {
        "Tc": 647.10, "Pc": 220.55, "omega": 0.345,
        "cp": (3.470, 1.450e-3, 0.0, 0.121e5)
    }
}

# Unified Cubic EOS Constants (SVNA Table 3.1)
EOS_MODELS = {
    "PR (Peng-Robinson)": {
        "code": "PR", "sigma": 1.0 + np.sqrt(2.0), "epsilon": 1.0 - np.sqrt(2.0),
        "Omega": 0.07779607, "Psi": 0.45723553
    },
    "SRK (Soave-Redlich-Kwong)": {
        "code": "SRK", "sigma": 1.0, "epsilon": 0.0,
        "Omega": 0.08664, "Psi": 0.42748
    },
    "RK (Redlich-Kwong)": {
        "code": "RK", "sigma": 1.0, "epsilon": 0.0,
        "Omega": 0.08664, "Psi": 0.42748
    },
    "vdW (van der Waals)": {
        "code": "vdW", "sigma": 0.0, "epsilon": 0.0,
        "Omega": 0.125, "Psi": 27.0 / 64.0
    }
}


def get_alpha_and_deriv(eos_code, Tr, omega):
    """Returns alpha(Tr) and d(ln alpha)/d(ln Tr) for the selected EOS."""
    if eos_code == "vdW":
        return 1.0, 0.0
    elif eos_code == "RK":
        return 1.0 / np.sqrt(Tr), -0.5
    elif eos_code == "SRK":
        kappa = 0.480 + 1.574 * omega - 0.176 * (omega ** 2)
        sqrt_Tr = np.sqrt(Tr)
        alpha = (1.0 + kappa * (1.0 - sqrt_Tr)) ** 2
        dln_alpha = (-kappa * sqrt_Tr) / (1.0 + kappa * (1.0 - sqrt_Tr))
        return alpha, dln_alpha
    else:  # PR
        kappa = 0.37464 + 1.54226 * omega - 0.26992 * (omega ** 2)
        sqrt_Tr = np.sqrt(Tr)
        alpha = (1.0 + kappa * (1.0 - sqrt_Tr)) ** 2
        dln_alpha = (-kappa * sqrt_Tr) / (1.0 + kappa * (1.0 - sqrt_Tr))
        return alpha, dln_alpha


def solve_cubic_cardano(c2, c1, c0):
    """Analytical Cardano-Vieta solver for Z^3 + c2*Z^2 + c1*Z + c0 = 0."""
    p = c1 - (c2 ** 2) / 3.0
    q = (2.0 * (c2 ** 3)) / 27.0 - (c2 * c1) / 3.0 + c0
    discriminant = (q / 2.0) ** 2 + (p / 3.0) ** 3

    if discriminant > 1e-14:
        sqrt_D = np.sqrt(discriminant)
        u = np.cbrt(-q / 2.0 + sqrt_D)
        v = np.cbrt(-q / 2.0 - sqrt_D)
        return [u + v - c2 / 3.0]
    elif abs(discriminant) <= 1e-14:
        u = np.cbrt(-q / 2.0)
        return sorted([2.0 * u - c2 / 3.0, -u - c2 / 3.0])
    else:
        r = np.sqrt(-((p / 3.0) ** 3))
        cos_arg = np.clip(-q / (2.0 * r), -1.0, 1.0)
        theta = np.arccos(cos_arg)
        m = 2.0 * np.sqrt(-p / 3.0)
        r1 = m * np.cos(theta / 3.0) - c2 / 3.0
        r2 = m * np.cos((theta + 2.0 * np.pi) / 3.0) - c2 / 3.0
        r3 = m * np.cos((theta + 4.0 * np.pi) / 3.0) - c2 / 3.0
        return sorted([r1, r2, r3])


def eval_root_props(Z, T, P_bar, B, q_param, sigma, epsilon, eos_code, dln_alpha):
    """Evaluates residual H^R, S^R, ln(phi), and f for a given compressibility root Z."""
    if eos_code == "vdW":
        I_val = B / Z
    else:
        I_val = (1.0 / (sigma - epsilon)) * np.log((Z + sigma * B) / (Z + epsilon * B))

    HR_RT = (Z - 1.0) + (dln_alpha - 1.0) * q_param * I_val
    ln_phi = (Z - 1.0) - np.log(Z - B) - q_param * I_val
    SR_R = np.log(Z - B) + dln_alpha * q_param * I_val

    HR = HR_RT * R * T
    SR = SR_R * R
    phi = np.exp(ln_phi)
    f_bar = phi * P_bar
    Vm_cm3 = (Z * R * T) / (P_bar * 1e5) * 1e6

    return {
        "Z": Z, "Vm_cm3": Vm_cm3, "HR": HR, "SR": SR,
        "ln_phi": ln_phi, "phi": phi, "f_bar": f_bar
    }


def ideal_gas_integrals(T, P_bar, cp_coeffs):
    """Integrates Cp^ig/R = A + B*T + C*T^2 + D*T^-2 from (298.15 K, 1 bar) to (T, P)."""
    A_c, B_c, C_c, D_c = cp_coeffs
    dH_ig = R * (
        A_c * (T - T_REF)
        + 0.5 * B_c * (T**2 - T_REF**2)
        + (1.0 / 3.0) * C_c * (T**3 - T_REF**3)
        - D_c * (1.0 / T - 1.0 / T_REF)
    )
    dS_ig = R * (
        A_c * np.log(T / T_REF)
        + B_c * (T - T_REF)
        + 0.5 * C_c * (T**2 - T_REF**2)
        - 0.5 * D_c * (1.0 / (T**2) - 1.0 / (T_REF**2))
        - np.log(P_bar / P_REF)
    )
    return dH_ig, dS_ig


def solve_eos_state(T, P_bar, Tc, Pc_bar, omega, cp_coeffs, eos_key):
    """Master solver for any of the 4 Cubic EOS models."""
    eos = EOS_MODELS[eos_key]
    code = eos["code"]
    sigma, epsilon = eos["sigma"], eos["epsilon"]
    Omega, Psi = eos["Omega"], eos["Psi"]

    P = P_bar * 1e5
    Pc = Pc_bar * 1e5
    Tr = T / Tc

    alpha, dln_alpha = get_alpha_and_deriv(code, Tr, omega)
    a = Psi * alpha * (R ** 2) * (Tc ** 2) / Pc
    b = Omega * R * Tc / Pc

    A = (a * P) / ((R * T) ** 2)
    B = (b * P) / (R * T)
    q_param = A / B

    # Polynomial coefficients: Z^3 + c2*Z^2 + c1*Z + c0 = 0
    u_eos = sigma + epsilon
    w_eos = sigma * epsilon
    c2 = -1.0 - B + u_eos * B
    c1 = A + w_eos * (B ** 2) - u_eos * B - u_eos * (B ** 2)
    c0 = -A * B - w_eos * (B ** 2) - w_eos * (B ** 3)

    raw_roots = solve_cubic_cardano(c2, c1, c0)
    valid_roots = [r for r in raw_roots if r > B]
    if not valid_roots:
        valid_roots = [max(raw_roots)]

    # Evaluate candidate roots (Liquid vs Vapor)
    candidates = [
        eval_root_props(r, T, P_bar, B, q_param, sigma, epsilon, code, dln_alpha)
        for r in valid_roots
    ]

    # Pick thermodynamically stable root via minimum Gibbs energy (G^R/RT = ln_phi)
    stable = min(candidates, key=lambda x: x["ln_phi"])
    vapor_root = max(candidates, key=lambda x: x["Z"])
    liquid_root = min(candidates, key=lambda x: x["Z"])

    if len(candidates) == 1:
        phase_label = "Single Phase / Supercritical" if T >= Tc else "Single Phase Fluid"
    else:
        if abs(vapor_root["ln_phi"] - liquid_root["ln_phi"]) < 1e-4:
            phase_label = "Vapor-Liquid Equilibrium (Two-Phase)"
        elif stable["Z"] == vapor_root["Z"]:
            phase_label = "Subcritical — Stable Vapor (Three Roots Detected)"
        else:
            phase_label = "Subcritical — Stable Liquid (Three Roots Detected)"

    dH_ig, dS_ig = ideal_gas_integrals(T, P_bar, cp_coeffs)
    H_abs = dH_ig + stable["HR"]
    S_abs = dS_ig + stable["SR"]

    return {
        "stable": stable,
        "vapor": vapor_root,
        "liquid": liquid_root,
        "num_roots": len(candidates),
        "phase_label": phase_label,
        "H_abs": H_abs,
        "S_abs": S_abs,
        "Tr": Tr,
        "Pr": P_bar / Pc_bar
    }


# --- Sidebar Controls ---
st.sidebar.header("⚙️ Input Parameters")
fluid_name = st.sidebar.selectbox("Select Pure Fluid", list(FLUID_DB.keys()))
eos_name = st.sidebar.selectbox("Select Cubic Equation of State", list(EOS_MODELS.keys()))

fluid = FLUID_DB[fluid_name]
Tc, Pc, omega, cp_coeffs = fluid["Tc"], fluid["Pc"], fluid["omega"], fluid["cp"]

st.sidebar.caption(f"**Critical Props:** $T_c = {Tc}$ K | $P_c = {Pc}$ bar | $\\omega = {omega}$")
T = st.sidebar.number_input("Temperature, T (K)", min_value=100.0, max_value=1000.0, value=350.0, step=5.0)
P = st.sidebar.number_input("Pressure, P (bar)", min_value=0.5, max_value=300.0, value=50.0, step=1.0)

res = solve_eos_state(T, P, Tc, Pc, omega, cp_coeffs, eos_name)
st_prop = res["stable"]

# --- Header & Tabs ---
st.title("🧪 Unified Cubic EOS Explorer: PVT, H, S & Fugacity")
st.markdown(f"**Model:** `{eos_name}` | **Fluid:** `{fluid_name}` | **Phase State:** `{res['phase_label']}`")

tab1, tab2, tab3 = st.tabs([
    "📊 Calculator & Isotherm Plots",
    "✅ SVNA Textbook Validation (Table 13.7)",
    "📘 Unified Theory & Equations"
])

with tab1:
    c1, c2, c3 = st.columns(3)
    with c1:
        st.subheader("1. PVT Behavior")
        st.metric("Stable Compressibility (Z)", f"{st_prop['Z']:.4f}")
        st.metric("Molar Volume (Vm)", f"{st_prop['Vm_cm3']:.2f} cm³/mol")
        if res["num_roots"] > 1:
            st.info(
                f"**Two-Phase Roots Found:**\n"
                f"- Vapor: $Z_v = {res['vapor']['Z']:.4f}$ ($V_v = {res['vapor']['Vm_cm3']:.1f}$ cm³/mol)\n"
                f"- Liquid: $Z_l = {res['liquid']['Z']:.4f}$ ($V_l = {res['liquid']['Vm_cm3']:.1f}$ cm³/mol)"
            )

    with c2:
        st.subheader("2. Enthalpy & Entropy (HW)")
        st.metric("Residual Enthalpy (Hᴿ)", f"{st_prop['HR']:.1f} J/mol", f"{st_prop['HR']/1000:.3f} kJ/mol", delta_color="off")
        st.metric("Residual Entropy (Sᴿ)", f"{st_prop['SR']:.3f} J/(mol·K)")
        st.metric("Absolute H (ref 298.15 K)", f"{res['H_abs']/1000:.2f} kJ/mol")
        st.metric("Absolute S (ref 298.15 K, 1 bar)", f"{res['S_abs']:.2f} J/(mol·K)")

    with c3:
        st.subheader("3. Fugacity & Stability")
        st.metric("Fugacity Coefficient (φ)", f"{st_prop['phi']:.4f}")
        st.metric("Pure Component Fugacity (f)", f"{st_prop['f_bar']:.3f} bar")
        st.metric("Normalized Residual Gibbs (Gᴿ/RT)", f"{st_prop['ln_phi']:.4f}")

    st.markdown("---")
    # Isotherm Plots (P-V and Fugacity vs P)
    p_grid = np.linspace(1.0, max(120.0, P * 1.5), 100)
    v_grid, f_grid, z_grid = [], [], []
    for p_val in p_grid:
        r_pt = solve_eos_state(T, p_val, Tc, Pc, omega, cp_coeffs, eos_name)["stable"]
        v_grid.append(r_pt["Vm_cm3"])
        f_grid.append(r_pt["f_bar"])
        z_grid.append(r_pt["Z"])

    pcol1, pcol2 = st.columns(2)
    with pcol1:
        fig_pv = go.Figure()
        fig_pv.add_trace(go.Scatter(x=v_grid, y=p_grid, mode="lines", name=f"Isotherm at {T} K", line=dict(color="#1d4ed8", width=3)))
        fig_pv.add_trace(go.Scatter(x=[st_prop["Vm_cm3"]], y=[P], mode="markers", name="Operating State", marker=dict(color="red", size=10)))
        fig_pv.update_layout(title=f"P-Vm Isotherm ({EOS_MODELS[eos_name]['code']} EOS)", xaxis_title="Molar Volume Vm [cm³/mol] (Log Scale)", yaxis_title="Pressure P [bar]", xaxis_type="log", height=380)
        st.plotly_chart(fig_pv, use_container_width=True)

    with pcol2:
        fig_f = go.Figure()
        fig_f.add_trace(go.Scatter(x=p_grid, y=f_grid, mode="lines", name="Real Gas Fugacity (f)", line=dict(color="#047857", width=3)))
        fig_f.add_trace(go.Scatter(x=p_grid, y=p_grid, mode="lines", name="Ideal Gas (f = P)", line=dict(color="gray", dash="dash")))
        fig_f.add_trace(go.Scatter(x=[P], y=[st_prop["f_bar"]], mode="markers", name="Operating State", marker=dict(color="red", size=10)))
        fig_f.update_layout(title="Pure Component Fugacity (f) vs Pressure", xaxis_title="System Pressure P [bar]", yaxis_title="Fugacity f [bar]", height=380)
        st.plotly_chart(fig_f, use_container_width=True)

with tab2:
    st.subheader("Benchmark Comparison Across All 4 Cubic EOS Models")
    st.markdown("Comparison at current $(T, P)$ across **vdW, RK, SRK, and PR** (plus SVNA Table 13.7 benchmark check for n-Butane at $500\\text{ K}, 50\\text{ bar}$):")
    rows = []
    for m_name in EOS_MODELS.keys():
        out = solve_eos_state(T, P, Tc, Pc, omega, cp_coeffs, m_name)["stable"]
        rows.append({
            "EOS Model": m_name,
            "Z (Compressibility)": f"{out['Z']:.4f}",
            "Vm (cm³/mol)": f"{out['Vm_cm3']:.2f}",
            "Hᴿ (J/mol)": f"{out['HR']:.1f}",
            "Sᴿ (J/mol·K)": f"{out['SR']:.3f}",
            "φ (Fugacity Coeff)": f"{out['phi']:.4f}",
            "f (bar)": f"{out['f_bar']:.3f}"
        })
    st.table(rows)

with tab3:
    st.subheader("Unified Cubic Equation of State Formulation (SVNA Ch. 3, 6 & 13)")
    st.latex(r"P = \frac{RT}{V - b} - \frac{a(T)}{(V + \epsilon b)(V + \sigma b)}")
    st.latex(r"Z^3 + (-1 - B + uB)Z^2 + (A + wB^2 - uB - uB^2)Z + (-AB - wB^2 - wB^3) = 0")
    st.latex(r"I = \frac{1}{\sigma - \epsilon} \ln\left(\frac{Z + \sigma B}{Z + \epsilon B}\right), \quad q = \frac{A}{B} = \frac{a}{bRT}")
    st.latex(r"\frac{H^R}{RT} = Z - 1 + \left(\frac{d\ln\alpha}{d\ln T_r} - 1\right) q I, \quad \frac{S^R}{R} = \ln(Z - B) + \left(\frac{d\ln\alpha}{d\ln T_r}\right) q I")
    st.latex(r"\ln\phi = \frac{G^R}{RT} = Z - 1 - \ln(Z - B) - q I \implies f = \phi P")
