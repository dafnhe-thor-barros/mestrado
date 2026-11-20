import os
import pypdf

bib_dir = r"C:\ZZ - Antigravity\Mestrado\projeto-mestrado\bibliografia"
files = [f for f in os.listdir(bib_dir) if f.endswith('.pdf')]

# Let's check which files exist and extract text from the recently downloaded ones to confirm their identity.
recent_files = [
    "Ali_et_al_2019_Towards_Smarter_Battery_Management_System_SOC.pdf",
    "Constant-Power_Load_System_Stabilization_by_Passive_Damping.pdf",
    "Current_collection_and_discharges_on_high_voltage_solar_array_for_space_use_in_plasma.pdf",
    "Development_of_400_V_Solar_Array_Technology_for_Low_Earth_Orbit_Plasma_Environment.pdf",
    "Study_of_the_Plasma_Interference_With_High-Voltage_Electrode_Array_for_Space_Power_Application.pdf",
    "A_Comparative_Study_on_Maximum_Power_Point_Tracking_Techniques_for_Photovoltaic_Power_Systems.pdf",
    "Impedance-Based_Stability_Criterion_for_Grid-Connected_Inverters.pdf",
    "Talebzadeh_Beik_2024_MVDC_Large_Spacecraft_Megawatt_Power_Propulsion.pdf",
    "A_primary_side_clamping_circuit_applied_to_the_ZVS-PWM_asymmetrical_half-bridge_converter.pdf"
]

print("Verifying recently downloaded files:")
for rf in recent_files:
    p = os.path.join(bib_dir, rf)
    if os.path.exists(p):
        try:
            reader = pypdf.PdfReader(p)
            text = reader.pages[0].extract_text()[:250].replace('\n', ' ')
            print(f"[OK] {rf} ({len(reader.pages)} pgs): {text[:120]}...")
        except Exception as e:
            print(f"[ERR] {rf}: {e}")
    else:
        print(f"[MISSING] {rf}")
