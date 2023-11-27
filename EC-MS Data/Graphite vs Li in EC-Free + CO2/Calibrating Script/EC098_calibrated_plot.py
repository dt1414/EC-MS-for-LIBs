import ixdat
from ixdat import Measurement
#from navigation import DATA_DIR
from pathlib import Path

# Ensure that all quantification uses the Spectro Inlets quantification package:
ixdat.config.plugins.USE_QUANT = True

#Data_Dir = Path("~/Users/daisythornton/Desktop/Quantification/EC021_Graphite_vs_Li")


# --- Import data --- #

# experiment data:
ec_ms = Measurement.read("/Users/daisythornton/Library/Mobile Documents/com~apple~CloudDocs/EH Backup/PhD/Data/EC-MS/2022/Papers/Transition Metal Measurments/EC098/my_combined_data.csv", reader="ixdat")
# calibration data:
ms_cal = Measurement.read(
    "/Users/daisythornton/Library/Mobile Documents/com~apple~CloudDocs/EH Backup/PhD/Data/EC-MS/2022/Papers/Methods Measurements/EC021 _Graphite_vs_Li/2022-06-15 10_38_59 Cal_Ar-bg_ethylene_mix.tsv",
    technique="MS",
)
ms = Measurement.read("/Users/daisythornton/Library/Mobile Documents/com~apple~CloudDocs/EH Backup/PhD/Data/EC-MS/2022/Papers/Transition Metal Measurments/EC098/2022-10-18 10_45_30 EC098.tsv",
    technique="MS",
)

we_area = 14.726 #working electrode area in mm2

# --- Plot the full MS data for Ar normalisation --- #

ms.plot() #grab tspan from this plot where Ar is stable and measured from glovebox only
S_M40_exp = ms.grab_signal("M40", tspan = [142, 1050])[1].mean()
S_M40_cal = ms_cal.grab_signal("M40", tspan = [400, 500])[1].mean()
factor = S_M40_exp / S_M40_cal #calculate factor by which to scale sensitivity factors
from ixdat.si_quant_patch import scale_by_factor

# --- Calculate the calibration --- #

# Direct addition (below) should work in the future. For now, though, use this:
from ixdat.quant_patch import append_sensitivity_factors
from spectro_inlets_quantification.medium import Medium

medium = Medium()
medium.p = 1.015e5   # pressure in Pa

# Use when matrix not working 
calibration = append_sensitivity_factors(
    ms_cal.gas_flux_calibration(mol="Ar", mass="M40", tspan=[15100, 15273]),
    ms_cal.gas_flux_calibration(mol="He", mass="M4", tspan=[16800, 16900]),
    ms_cal.multicomp_gas_flux_calibration(
        gas={"Ar": 1, "CO2": 1140e-6, "O2": 1053e-6, "H2": 931e-6, "C2H4": 877e-6, "CO": 991e-6},
        tspan=[18400, 18524],
        gas_bg="Ar",
        tspan_bg=[15100, 15273],
        mol_list=[
            "CO2", "O2",  "H2", "C2H4",
            "CO"
        ],
        mass_list=[
            "M44", "M32",  "M2", "M26",
            "M28"
        ],
    ),
)

#ms_cal.plot()

calibration = scale_by_factor(calibration, factor=factor)


#Predicting solvent contribution - good to check if significant difference 
#fit = calibration.fit

#fit.plot_F_vs_f(
#    predict={"CO": "M28", "ethylene_carbonate": "M29", "ethyl_methyl_carbonate": "M45"}
#)

#sf_EC_M29 = calibration.fit.predict_sf(mol="ethylene_carbonate", mass="M29")
#sf_EMC_M45 = calibration.fit.predict_sf(mol="ethyl_methyl_carbonate", mass="M45")

#calibration.append(sf_EC_M29)
#calibration.append(sf_EMC_M45)

#calibration = (
#     ms_cal.gas_flux_calibration(mol="Ar", mass="M40", tspan=[15100, 15273])
#    + ms_cal.gas_flux_calibration(mol="He", mass="M4", tspan=[16800, 16900])
#    + ms_cal.multicomp_gas_flux_calibration(
#        gas={"Ar": 1, "CO2": 1140e-6, "O2": 1053e-6, "H2": 931e-6, "C2H4": 877e-6, "CO": 991e-6},
#        tspan=[18400, 18524],
#        gas_bg="Ar",
#        tspan_bg=[15100, 15273],
#        mol_list=["CO2", "O2", "H2", "C2H4", "CO"],
#        mass_list=["M44", "M32", "M2", "M26", "M28"],
#    )
# )

# --- Adding isotopes - This code will be replaced with a quant update --- #
#from ixdat.quant_patch import add_isotopes

#add_isotopes(
#    calibration, {"CO2": ("M44", ["M46", "M48"]), "O2": ("M32", ["M34", "M36"])}
#)
# # ^^---vv After a quant update the above code should be replaced by this: ^^---vv
# calibration.add_isotopes(
#    {"CO2": ("M44", ["M46", "M48"]), "O2": ("M32", ["M34", "M36"])}
# )
# # ---

#calibration.save("My isotope calibration")

ec_ms.set_quantifier(calibration=calibration, carrier="He")

# --- Plot calibrated data with ixdat's plotting function --- #

ec_ms.plot(
    mol_list=["O2", "CO2", "H2", "C2H4", "CO"],
    logplot=False,
    #tspan_bg=[135000, 135000],
    unit="pmol/s"
)

ec_ms.plot(
    mol_list=["H2"],
    logplot=False,
    unit="pmol/s"
    )

ec_ms.plot(
    mol_list=["C2H4"],
    logplot=False,
    unit="pmol/s"
    )

ec_ms.plot(
    mol_list=["CO"],
    logplot=False,
    unit="pmol/s"
    )

ec_ms.plot(
    mol_list=["CO2"],
    logplot=False,
    unit="pmol/s"
    )

# - Manually plot calibrated data on three axes with smoothing and moving background - #


# The code for setting up the axes is adapted from here:
#   https://github.com/ScottSoren/Huang2021/blob/main/Figure%203/fig3.py
from matplotlib import pyplot as plt, gridspec
from ixdat.plotters.plotting_tools import smooth_vector, calc_linear_background, color_axis

# prepare the axes using a gridspec:
fig_a = plt.figure()

gs_a = gridspec.GridSpec(4, 1, fig_a)
# gs.update(hspace=0.025)
ax_a_H2 = plt.subplot(gs_a[0, 0])  # H2 flux axis
ax_a_CO2 = plt.subplot(gs_a[1, 0])  # CO flux axis
ax_a_C2H4 = plt.subplot(gs_a[2, 0]) # C2H4 flux axis
ax_a_U = plt.subplot(gs_a[3, 0])  # potential axis
ax_a_J = ax_a_U.twinx()  # current axis

fig_a.set_figheight(fig_a.get_figwidth() * 1.25)

tspan = [0, 45000]
tspans_bg_H2 = [[79, 80], [463, 464], [654, 655], [878, 879], [974, 975], [1133, 1134], [1293, 1294], [1517, 1518], [1708, 1709], [1806, 1807], [1900, 1901], [2187, 2188], [2357, 2357], [2411, 2412], [2629, 2630], [2730, 2731], [2857, 2858], [3018, 3019], [3093, 3094], [3347, 3348], [3401, 3402], [3575, 3576], [3752, 3753], [4104, 4105], [6052, 6053], [6531, 6532], [6851, 6852], [7425, 7426], [7681, 7682], [8352, 8353], [8927, 8928], [9565, 9566], [10332, 10333], [10747, 10748], [11929, 11930], [12504, 12505], [13175, 13176], [13845, 13846], [14324, 14325]
                #[30403, 30404], [30535, 30536], [30668, 30669], [30800, 30801], [33502, 33503], [35171, 35172], [40840, 40841]
                ]
tspans_bg_C2H4 = [[79, 80], [463, 464], [654, 655], [878, 879], [974, 975], [1133, 1134], [1293, 1294], [1517, 1518], [1708, 1709], [1806, 1807], [1900, 1901], [2187, 2188], [2357, 2357], [2411, 2412], [2629, 2630], [2730, 2731], [2857, 2858], [3018, 3019], [3093, 3094], [3347, 3348], [3401, 3402], [3575, 3576], [3752, 3753], [4104, 4105], [6052, 6053], [6531, 6532], [6851, 6852], [7425, 7426], [7681, 7682], [8352, 8353], [8927, 8928], [9565, 9566], [10332, 10333], [10747, 10748], [11929, 11930], [12504, 12505], [13175, 13176], [13845, 13846], [14324, 14325]
                  #[4961, 4962], [8039, 8040], [10541, 10542], [13491, 13492], [16634, 16635], [19520, 19521], [22599, 22600], [25614, 25615], [28692, 28693], [31386, 31387], [36774, 36775]
                  ]
#tspans_bg_CO = [[29, 30], [42, 43], [460, 461], [878, 879], [1244, 1245], [1558, 1559], [2029, 2030], [2866, 2867], [3098, 3099], [4060, 4061], [4472, 4473], [7952, 7953], [8960, 8961], [11021, 11022], [11937, 11938], [13998, 13999], [15005, 15006], [16929, 16930], [17891, 17892], [19906, 19907], [20913, 20914], [22791, 22792], [23982, 23983], [31768, 31769], [33004, 33005], [34790, 34791], [35615, 35616], [36485, 36486], [37218, 37219], [37630, 37631], [38317, 38318], [39462, 39463], [40653, 40654]
                #, [4042, 4043], [8342, 8343], [11577, 11578], [14615, 14616], [17613, 17614], [20611, 20612], [23728, 23729], [26686, 26687], [29882, 29883], [30355, 30356], [31184, 31185], [31341, 31342], [32762, 32763], [35641, 35642], [37772, 37773], [40730, 40731]
#                ]

# Plot the potential and current:
t_U, U = ec_ms.grab("potential", tspan=tspan)
t_J, J = ec_ms.grab("current", tspan=tspan)
ax_a_U.plot(t_U, U, color="black")
ax_a_J.plot(t_J, J / we_area, color="red")
ax_a_U.set_ylabel("Potential (V vs. Li)")
ax_a_J.set_ylabel("Current Density (mA mm$^{-2}$)")
color_axis(ax_a_J, "red")
ax_a_U.set_xlabel("Time (s)")
ax_a_U.set_xlim(t_U[0], t_U[-1])
ax_a_J.set_xlim(t_U[0], t_U[-1])

t, n_dot_dict = ec_ms.grab_fluxes(tspan=tspan)

# for named colors, see: https://matplotlib.org/stable/gallery/color/named_colors.html

background_subtracted_fluxes = {}

# Plot the H2:
for mol, color in [("H2", "royalblue")]:
    n_dot = n_dot_dict[mol]
    n_dot = smooth_vector(n_dot, 10)
    n_dot_bg = calc_linear_background(t, n_dot, tspans=tspans_bg_H2)
    ax_a_H2.plot(t, (n_dot - n_dot_bg) * (1e12 / we_area), color=color)
    background_subtracted_fluxes[mol] = (n_dot - n_dot_bg) * (1e12 / we_area)
ax_a_H2.tick_params()
ax_a_H2.set_ylabel("Gas Evolution Rate (pmol s$^{-1}$ mm$^{-2}$)                                ")
ax_a_H2.set_xlabel("Time (s)")
ax_a_H2.xaxis.set_label_position("top")
ax_a_H2.tick_params(
    axis="x", top=True, bottom=False, labeltop=True, labelbottom=False
)
ax_a_H2.set_xlim(t_U[0], t_U[-1])
ax_a_H2.set_ylim(-0.05, 0.5)

# Plot the CO2:
for mol, color in [("CO2", "firebrick")]:
    n_dot = n_dot_dict[mol]
    n_dot = smooth_vector(n_dot, 10)
    #n_dot_bg = calc_linear_background(t, n_dot, tspans=tspans_bg_CO)
    ax_a_CO2.plot(t, (n_dot) * (1e12 / we_area), color=color)
    background_subtracted_fluxes[mol] = (n_dot) * (1e12 / we_area)
ax_a_CO2.tick_params(axis="x", top=True, bottom=True, labeltop=False, labelbottom=False)
ax_a_CO2.set_xlabel("")
ax_a_CO2.set_xlim(t_U[0], t_U[-1])
ax_a_CO2.set_yscale("log")
ax_a_CO2.set_ylim(0.6, 400)

# Plot C2H4:
for mol, color in [("C2H4", "darkgreen")]:
    n_dot = n_dot_dict[mol]
    n_dot = smooth_vector(n_dot, 10)
    n_dot_bg = calc_linear_background(t, n_dot, tspans=tspans_bg_C2H4)
    ax_a_C2H4.plot(t, (n_dot - n_dot_bg) * (1e12 / we_area), color=color)
    background_subtracted_fluxes[mol] = (n_dot - n_dot_bg) * (1e12 / we_area)
ax_a_C2H4.tick_params(axis="x", top=True, bottom=True, labeltop=False, labelbottom=False)
ax_a_C2H4.set_xlabel("")
ax_a_C2H4.set_xlim(t_U[0], t_U[-1])
ax_a_C2H4.set_ylim(-0.25, 3)

background_subtracted_fluxes["t"] = t
background_subtracted_fluxes["U"] = ec_ms.grab_for_t("potential", t)
background_subtracted_fluxes["J"] = ec_ms.grab_for_t("current", t)

fig_a.savefig("/Users/daisythornton/Library/Mobile Documents/com~apple~CloudDocs/EH Backup/PhD/Data/EC-MS/2022/Papers/Transition Metal Measurments/EC098/EC098.svg")

from pandas import DataFrame

#df = DataFrame(background_subtracted_fluxes)
#df.to_csv(DATA_DIR / "results.csv")


