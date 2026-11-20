import os
import pypdf

bib_dir = r"C:\ZZ - Antigravity\Mestrado\projeto-mestrado\bibliografia"
files = [f for f in os.listdir(bib_dir) if f.endswith('.pdf')]

# Let's map each of the 26 references
refs = [
    (1, "ALI, M. U. et al. (2019)", "Towards a smarter battery management system for electric vehicle applications: A critical review of lithium-ion battery state of charge estimation", "Ali_et_al_2019_Towards_Smarter_Battery_Management_System_SOC.pdf", "https://doi.org/10.3390/en12030446"),
    (2, "CESPEDES, M.; XING, L.; SUN, J. (2011)", "Constant-power load system stabilization by passive damping", "Cespedes_Xing_Sun_2011_Constant_Power_Load_Stabilization_Passive_Damping.pdf", "https://doi.org/10.1109/TPEL.2011.2151880"),
    (3, "CHEN, Haonan et al. (2026)", "Edge intelligence for space-air-ground integrated network", "Chen_2026_Edge_Intelligence_Space_Air_Ground_Network.pdf", "https://doi.org/10.1016/j.cja.2025.103969"),
    (4, "FORTESCUE, Peter; SWINERD, Graham; STARK, John (2011)", "Spacecraft systems engineering (4. ed.)", "Fortescue_Swinerd_Stark_2011_Spacecraft_Systems_Engineering.pdf", "https://www.wiley.com/en-us/Spacecraft+Systems+Engineering%2C+4th+Edition-p-9780470750124"),
    (5, "FUJII, H.; ABE, T. (2003)", "Current collection and discharges on high voltage solar array for space use in plasma", "Fujii_Abe_2003_Current_Collection_Discharges_High_Voltage_Solar_Array.pdf", "https://doi.org/10.1109/TPS.2003.818448"),
    (6, "HILLARD, G. B.; FERGUSON, D. C. (2003)", "LEO Spacecraft Charging Design Guidelines (NASA/TM-2003-212737)", "Hillard_Ferguson_2003_LEO_Spacecraft_Charging_Guidelines_NASA_TM.pdf", "https://ntrs.nasa.gov/citations/20040026971"),
    (7, "HOSODA, S. et al. (2006)", "Development of 400 V solar array technology for low earth orbit plasma environment", "Hosoda_et_al_2006_Development_400V_Solar_Array_LEO_Plasma.pdf", "https://doi.org/10.1109/TPS.2006.883287"),
    (8, "IWASA, Minoru et al. (2006)", "Study of the plasma interference with high-voltage electrode array for space power application", "Iwasa_et_al_2006_Plasma_Interference_High_Voltage_Electrode_Array.pdf", "https://doi.org/10.1109/TPS.2006.883369"),
    (9, "KOPP, G.; LEAN, J. L. (2011)", "A new, lower value of total solar irradiance", "Kopp_Lean_2011_New_Lower_Value_Total_Solar_Irradiance_GRL.pdf", "https://doi.org/10.1029/2010GL045777"),
    (10, "LYDEN, Sarah; HAQUE, Md. Enamul (2015)", "Maximum Power Point Tracking techniques for photovoltaic systems: A review", "Lyden_Haque_2015_Maximum_Power_Point_Tracking_Techniques_Review.pdf", "https://doi.org/10.1016/j.rser.2015.07.172"),
    (11, "MALIK, M. Z. et al. (2020)", "Investigation of MPPT techniques under uniform and non-uniform solar irradiation condition", "Malik_et_al_2020_Investigation_MPPT_Uniform_NonUniform_Irradiation.pdf", "https://doi.org/10.1109/ACCESS.2020.3007710"),
    (12, "MIDDLEBROOK, R. D. (1976)", "Input filter considerations in design and application of switching regulators", "Middlebrook_1976_Input_Filter_Considerations_Switching_Regulators_PESC.pdf", "https://ieeexplore.ieee.org/document/5242354"),
    (13, "ORTIZ DEL CASTILLO, M. et al. (2024)", "Unlocking the potential of small satellites: Onboard AI processing on the SpIRIT CubeSat", "OrtizDelCastillo_2024_SpIRIT_Jetson.pdf", "https://arxiv.org/abs/2407.14031"),
    (14, "PATEL, M. R. (2004)", "Spacecraft power systems (CRC Press)", "Patel_2004_Spacecraft_Power_Systems.pdf", "https://www.routledge.com/Spacecraft-Power-Systems/Patel/p/book/9780849320866"),
    (15, "PLETT, G. L. (2004)", "Extended Kalman filtering for battery management systems of LiPB-based HEV battery packs: Part 1. Background", "Plett_2004_Extended_Kalman_Filtering_Battery_Management_HEV_LiPB.pdf", "https://doi.org/10.1016/j.jpowsour.2004.02.031"),
    (16, "RADWAN, A. A. A.; MOHAMED, Y. A.-R. I. (2011)", "Linear active stabilization of converter-dominated DC microgrids", "Radwan_Mohamed_2011_Linear_Active_Stabilization_DC_Microgrids.pdf", "https://doi.org/10.1109/TSG.2010.2093035"),
    (17, "RADWAN, A. A. A.; MOHAMED, Y. A.-R. I. (2012)", "Modeling, analysis, and stabilization of converter-fed AC microgrids", "Radwan_Mohamed_2012_Modeling_Analysis_Stabilization_AC_Microgrids.pdf", "https://doi.org/10.1109/TSG.2012.2191986"),
    (18, "RICCOBONO, A.; SANTI, E. (2014)", "Comprehensive review of stability criteria for DC power distribution systems with constant power loads", "Riccobono_Santi_2014_Comprehensive_Review_Stability_Criteria_DC_CPL.pdf", "https://doi.org/10.1109/TIA.2014.2309881"),
    (19, "SUBUDHI, B.; PRADHAN, R. (2013)", "A comparative study on maximum power point tracking techniques for photovoltaic power systems", "Subudhi_Pradhan_2013_Comparative_Study_MPPT_Photovoltaic.pdf", "https://doi.org/10.1109/TSTE.2012.2202294"),
    (20, "SUN, J. (2011)", "Impedance-based stability criterion for grid-connected inverters", "Sun_2011_Impedance_Based_Stability_Criterion_Grid_Connected_Inverters.pdf", "https://doi.org/10.1109/TPEL.2011.2136439"),
    (21, "TALEBZADEH, Sarah; BEIK, Omid (2024)", "Medium voltage DC system for large spacecraft megawatt-scale power and propulsion systems", "Talebzadeh_Beik_2024_MVDC_Large_Spacecraft_Megawatt_Power_Propulsion.pdf", "https://doi.org/10.3390/electronics13101810"),
    (22, "THALES ALENIA SPACE (2024)", "ASCEND: Advanced Space Cloud for European Net zero Emission and Data Sovereignty", "Thales_Alenia_Space_2024_ASCEND_Feasibility_Study_Space_Data_Centers.pdf", "https://www.thalesaleniaspace.com/en/press-release/ascend-study-shows-feasibility-space-data-centers-and-their-contribution-european-green"),
    (23, "VAN DER MERWE, J. et al. (2001)", "High voltage solar array plasma interactions: LEO and GEO environment simulations", "Van_der_Merwe_2001_High_Voltage_Solar_Array_Plasma_Interactions_SCTC.pdf", "https://ui.adsabs.harvard.edu/abs/2001sctc.conf..283V/abstract"),
    (24, "VARSHNI, Y. P. (1967)", "Temperature dependence of the energy gap in semiconductors", "Varshni_1967_Temperature_Dependence_Energy_Gap_Semiconductors.pdf", "https://doi.org/10.1016/0031-8914(67)90062-6"),
    (25, "WAN, H. et al. (2000)", "Active stabilization of constant power loads in DC distributed power systems", "Wan_et_al_2000_Active_Stabilization_Constant_Power_Loads_DC_DPS_PESC.pdf", "https://ieeexplore.ieee.org/document/878848"),
    (26, "XIONG, R. et al. (2018)", "Towards a smarter battery management system: A critical review on battery state of health monitoring methods", "Xiong_et_al_2018_Towards_Smarter_BMS_Critical_Review_Battery_SOH.pdf", "https://doi.org/10.1016/j.jpowsour.2018.10.019")
]

present_count = 0
pending_count = 0

print("=== STATUS DAS 26 REFERÊNCIAS DO ARTIGO 02 ===")
for num, author, title, expected_fn, link in refs:
    exists = os.path.exists(os.path.join(bib_dir, expected_fn))
    if exists:
        present_count += 1
        print(f"[{num:02d}] [PRESENTE] {author}")
        print(f"     Arquivo: {expected_fn}")
    else:
        pending_count += 1
        print(f"[{num:02d}] [PENDENTE] {author}")
        print(f"     Link: {link}")
        print(f"     Nome sugerido: {expected_fn}")

print(f"\nResumo: {present_count} PRESENTES, {pending_count} PENDENTES.")
