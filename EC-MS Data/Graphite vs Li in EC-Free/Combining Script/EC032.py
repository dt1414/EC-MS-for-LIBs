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
ms_data_dir = Path(r"D:\PhD\Data\EC-MS\2022\Papers\Methods Measurements\EC032 - Graphite vs Li (EC Free)")
ec_data_dir = Path(r"D:\PhD\Data\EC-MS\2022\Papers\Methods Measurements\EC032 - Graphite vs Li (EC Free)\EC032")

#Define mass spec file
ms = Measurement.read(ms_data_dir / "2022-05-12 17_18_23 EC_free_graphite.tsv", reader="zilien", technique="MS")
#ms = Measurement.read(ms_data_dir / "tmp", reader="zilien_tmp", technique="MS")
#ms.plot()

#Define electrochemistry file
ec = Measurement.read_set(ec_data_dir / "EC032_EC_free_G_vs_Li_02_CV_C01", suffix=".mpt", reader="biologic")

ec_ms = ec + ms

#Plot or export combined data
#ec.plot()
#ms.plot()
#ec_ms.plot()
ec_ms.plot(mass_list=["M2", "M26", "M27", "M28"])

ec_ms.export("EC032_combined_data.csv")