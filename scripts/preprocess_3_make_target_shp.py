"""
Author: Wenyu Ouyang
Date: 2025-02-03 06:34:39
LastEditTime: 2025-02-03 08:11:03
LastEditors: Wenyu Ouyang
Description: make shapefile for target stations
FilePath: \hydrotopo\scripts\preprocess_3_make_target_shp.py
Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
"""

import os
import geopandas as gpd
from pathlib import Path


project_dir = Path(os.path.abspath(__file__)).parent.parent
data_dir = os.path.join(project_dir, "data")
result_dir = project_dir.joinpath("results")
node_file = os.path.join(data_dir, "rr_stations", "rr_stations.shp")
nodes_gpd = gpd.read_file(node_file)
# 9 main reservoirs in Liaoning: 石佛寺，柴河，清河，闹德海，大伙房，观音阁，葠窝水库，汤河水库，白石水库
target_stcd_lst = [
    "20600340",
    "20800900",
    "20810200",
    "20910930",
    "21100150",
    "21110150",
    "21110400",
    "21113800",
    "21200510",
]
# find the row index of the target node with its STCD
target_indices = nodes_gpd[nodes_gpd["STCD"].isin(target_stcd_lst)].index
print(target_indices)
# make a new shape file with the target nodes
target_nodes_gpd = nodes_gpd.loc[target_indices]
target_dir = result_dir.joinpath("target_stations")
if not target_dir.exists():
    target_nodes_gpd.to_file(target_dir, encoding="utf-8")
