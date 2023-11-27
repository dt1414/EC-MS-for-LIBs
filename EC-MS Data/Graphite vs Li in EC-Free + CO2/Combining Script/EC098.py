# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:56:02 2022

@author: Goodenough
"""


from pathlib import Path
from ixdat import Measurement
from ixdat.techniques import MSMeasurement
from ixdat.techniques import ECMSMeasurement

#Define directory (folder containing data)
ms_data_dir = Path("/Users/daisythornton/Library/Mobile Documents/com~apple~CloudDocs/EH Backup/PhD/Data/EC-MS/2022/Papers/Transition Metal Measurments/EC098")
ec_data_dir = Path("/Users/daisythornton/Library/Mobile Documents/com~apple~CloudDocs/EH Backup/PhD/Data/EC-MS/2022/Papers/Transition Metal Measurments/EC098")

#Define mass spec file
ms = Measurement.read(ms_data_dir / "2022-10-18 10_45_30 EC098.tsv", reader="zilien", technique="MS")
#ms = Measurement.read(ms_data_dir / "tmp", reader="zilien_tmp", technique="MS")
#ms.plot()

#Define electrochemistry file
ec1 = Measurement.read_set(ec_data_dir / "EC098_EC_Free_CO2_01_OCV_C01", suffix=".mpt", reader="biologic")
ec2 = Measurement.read_set(ec_data_dir / "EC098_EC_Free_CO2_02_CV_C01", suffix=".mpt", reader="biologic")
ec3 = Measurement.read_set(ec_data_dir / "EC098_EC_Free_CO2_03_OCV_C01", suffix=".mpt", reader="biologic")
#ec4 = Measurement.read_set(ec_data_dir / "EC068_2 Cu spiked G vs Li_01_OCV_C01", suffix=".mpt", reader="biologic")
#ec5 = Measurement.read_set(ec_data_dir / "EC068_2 Cu spiked G vs Li_02_CV_C01", suffix=".mpt", reader="biologic")
#ec6 = Measurement.read_set(ec_data_dir / "EC068_2 Cu spiked G vs Li_03_OCV_C01", suffix=".mpt", reader="biologic")

ec_ms = (ec1 + ec2 + ec3) + ms

#Plot or export combined data
#ec.plot()
#ms.plot()
#ec_ms.plot(tspan=[0, 40000], mass_list=["M2", "M26", "M27", "M28", "M29", "M30"])
ec_ms.plot(mass_list=["M2", "M26", "M27", "M28", "M44", "M15"])

#ec_ms.export("my_combined_data.csv")
#ms.export("ms_data.csv")
#ec_ms.savefig("/Volumes/T7/PhD/Data/EC-MS/2022/Papers/Plots/EC080.svg")

