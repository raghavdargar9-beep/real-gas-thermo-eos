import streamlit as st
import numpy as np
import plotly.graph_objects as go

# Page Configuration
st.set_page_config(
    page_title="Real Gas Thermodynamics | PR-EOS",
    page_icon="🧪",
    layout="wide"
)

# Universal Gas Constant (J / mol·K)
R = 8.314462618

# Preset Gas Database: {Name: (Tc [K], Pc [bar], omega)}
GAS_DB = {
    "Carbon Dioxide (CO₂)": (304.12, 73.74, 0.224),
    "Methane (CH₄)": (190.56, 45.99, 0.011),
    "Nitrogen (N₂)": (126.20, 33.98, 0.037),
    "Ethane (C₂H₆)": (305.32, 48.72, 0.099),
    "Propane (C₃H₈)": (369.83, 42.48, 0.152),
    "Water / Steam (H₂O)": (647.14, 220.64, 0.344),
    "Custom Gas": (300.0, 50.0, 0.100)
}

def calculate_pr_eos(T, P_bar, Tc, Pc_bar, omega):
    """Solves Peng-Robinson EOS for PVT, H_dep, S_dep, and Fugacity."""
    P = P_bar * 1e5
    Pc = Pc_bar * 1e5

    # PR Constants
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * (omega ** 2)
    ac = 0.45723553 * (R ** 2) * (Tc ** 2) / Pc
    b = 0.07779607 * R * Tc / Pc

    Tr = T / Tc
    alpha = (1 + kappa * (1 - np.sqrt(Tr))) ** 2
    a = ac * alpha
    da_dT = -ac * kappa * np.sqrt(alpha / (T * Tc))

    A = (a * P) / ((R ** 2) * (T ** 2))
    B = (b * P) / (R * T)

    # Cubic Equation: Z^3 + c2*Z^2 + c1*Z + c0 = 0
    c2 = -(1.0 - B)
    c1 = A - 2.0 * B - 3.0 * (B ** 2)
    c0 = -(A * B - (B ** 2) - (B ** 3))

    roots = np.roots([1.0, c2, c1, c0])
    real_roots = roots[np.isreal(roots)].real
    # Filter positive roots greater than B
    valid_roots = [r for r in real_roots if r > B]
    Z = max(valid_roots) if valid_roots else max(real_roots)

    # Molar Volume (m^3/mol and L/mol)
    V_m3 = (Z * R * T) / P
    V_L = V_m3 * 1000.0
    V_ig_L = ((R * T) / P) * 1000.0

    # Logarithmic term for Departure & Fugacity
    sqrt2 = np.sqrt(2.0)
    log_term = np.log((Z + (1.0 + sqrt2) * B) / (Z + (1.0 - sqrt2) * B))

    # Homework: Enthalpy & Entropy Departure
    H_dep = R * T * (Z - 1.0) + ((T * da_dT - a) / (2.0 * sqrt2 * b)) * log_term
    S_dep = R * np.log(Z - B) + (da_dT / (2.0 * sqrt2 * b)) * log_term

    # Computer Project: Pure Component Fugacity
    ln_phi = (Z - 1.0) - np.log(Z - B) - (A / (2.0 * sqrt2 * B)) * log_term
    phi = np.exp(ln_phi)
    f_bar = phi * P_bar

    return {
        "Z": Z, "V_L": V_L, "V_ig_L": V_ig_L,
        "H_dep": H_dep, "S_dep": S_dep,
        "phi": phi, "f_bar": f_bar,
        "A": A, "B": B
    }

# --- Sidebar Inputs ---
st.sidebar.header("⚙️ Input Parameters")
selected_gas = st.sidebar.selectbox("Select Pure Component", list(GAS_DB.keys()))
default_Tc, default_Pc, default_omega = GAS_DB[selected_gas]

is_custom = (selected_gas == "Custom Gas")
Tc = st.sidebar.number_input("Critical Temp, Tc (K)", value=default_Tc, disabled=not is_custom, format="%.2f")
Pc = st.sidebar.number_input("Critical Pressure, Pc (bar)", value=default_Pc, disabled=not is_custom, format="%.2f")
omega = st.sidebar.number_input("Acentric Factor (ω)", value=default_omega, disabled=not is_custom, format="%.4f")

st.sidebar.markdown("---")
st.sidebar.subheader("Operating Conditions")
T = st.sidebar.slider("Temperature, T (K)", min_value=100.0, max_value=800.0, value=350.0, step=5.0)
P = st.sidebar.slider("Pressure, P (bar)", min_value=1.0, max_value=250.0, value=50.0, step=1.0)

# Run Calculation
res = calculate_pr_eos(T, P, Tc, Pc, omega)

# --- Main Dashboard ---
st.title("🧪 Real Gas Thermodynamics Simulator")
st.markdown("**Peng-Robinson Cubic Equation of State** | PVT Behavior, Departure Properties ($H, S$), and Pure Component Fugacity ($f$)")

# 3-Column Metric Layout for the 3 Assignment Parts
col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("1. PVT Behavior")
    st.metric("Compressibility Factor (Z)", f"{res['Z']:.5f}", delta=f"{res['Z']-1:.4f} vs Ideal")
    st.metric("Real Molar Volume (V)", f"{res['V_L']:.4f} L/mol")
    st.metric("Ideal Molar Volume (V_ig)", f"{res['V_ig_L']:.4f} L/mol")

with col2:
    st.subheader("2. Enthalpy & Entropy (HW)")
    st.metric("Enthalpy Departure (H - H_ig)", f"{res['H_dep']/1000:.3f} kJ/mol", f"{res['H_dep']:.1f} J/mol", delta_color="off")
    st.metric("Entropy Departure (S - S_ig)", f"{res['S_dep']:.4f} J/(mol·K)")
    st.caption("Residual properties relative to Ideal Gas state at (T, P).")

with col3:
    st.subheader("3. Pure Fugacity (Project)")
    st.metric("Fugacity Coefficient (φ)", f"{res['phi']:.5f}")
    st.metric("Fugacity (f)", f"{res['f_bar']:.3f} bar", delta=f"{res['f_bar']-P:.2f} bar vs P")
    st.metric("Reduced State (Tr, Pr)", f"Tr = {T/Tc:.2f}, Pr = {P/Pc:.2f}")

st.markdown("---")

# --- Interactive Isotherm Plots ---
st.subheader(f"📈 Isotherm Behavior at T = {T} K ({selected_gas})")
p_range = np.linspace(1.0, 250.0, 100)
z_vals, phi_vals, h_vals = [], [], []

for p_val in p_range:
    r = calculate_pr_eos(T, p_val, Tc, Pc, omega)
    z_vals.append(r["Z"])
    phi_vals.append(r["phi"])
    h_vals.append(r["H_dep"] / 1000.0)

plot_col1, plot_col2 = st.columns(2)

with plot_col1:
    fig_z = go.Figure()
    fig_z.add_trace(go.Scatter(x=p_range, y=z_vals, mode="lines", name="Z (PR-EOS)", line=dict(color="#2563eb", width=3)))
    fig_z.add_trace(go.Scatter(x=[P], y=[res["Z"]], mode="markers", name="Operating Point", marker=dict(color="red", size=10)))
    fig_z.add_hline(y=1.0, line_dash="dash", line_color="gray", annotation_text="Ideal Gas (Z=1)")
    fig_z.update_layout(title="Compressibility Factor (Z) vs Pressure", xaxis_title="Pressure (bar)", yaxis_title="Z", height=380)
    st.plotly_chart(fig_z, use_container_width=True)

with plot_col2:
    fig_f = go.Figure()
    fig_f.add_trace(go.Scatter(x=p_range, y=phi_vals, mode="lines", name="φ (Fugacity Coeff)", line=dict(color="#059669", width=3)))
    fig_f.add_trace(go.Scatter(x=[P], y=[res["phi"]], mode="markers", name="Operating Point", marker=dict(color="red", size=10)))
    fig_f.add_hline(y=1.0, line_dash="dash", line_color="gray", annotation_text="Ideal Gas (φ=1)")
    fig_f.update_layout(title="Fugacity Coefficient (φ) vs Pressure", xaxis_title="Pressure (bar)", yaxis_title="Fugacity Coefficient (φ)", height=380)
    st.plotly_chart(fig_f, use_container_width=True)

# --- Mathematical Formulation Expander ---
with st.expander("📘 View Governing Peng-Robinson EOS Equations & Derivations"):
    st.latex(r"P = \frac{RT}{V - b} - \frac{a(T)}{V(V + b) + b(V - b)}")
    st.latex(r"Z^3 - (1 - B)Z^2 + (A - 2B - 3B^2)Z - (AB - B^2 - B^3) = 0")
    st.latex(r"\frac{H - H^{ig}}{RT} = (Z - 1) + \frac{T\left(\frac{da}{dT}\right) - a}{2\sqrt{2}bRT} \ln\left[\frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B}\right]")
    st.latex(r"\frac{S - S^{ig}}{R} = \ln(Z - B) + \frac{\frac{da}{dT}}{2\sqrt{2}bR} \ln\left[\frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B}\right]")
    st.latex(r"\ln \phi = (Z - 1) - \ln(Z - B) - \frac{A}{2\sqrt{2}B} \ln\left[\frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B}\right]")
