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
# if we need to update the location of some points, we can use the following code
node_file_ndh = os.path.join(data_dir, "ndh", "ndh.shp")
if os.path.exists(node_file_ndh):
    nodes_gpd_ndh = gpd.read_file(node_file_ndh)
    # 明确列名
    source_id_field = 'id'
    target_id_field = 'STCD'
    # 确保ID字段的类型一致，转换为字符串类型
    nodes_gpd_ndh[source_id_field] = nodes_gpd_ndh[source_id_field].astype(str)
    # 创建一个字典，键为ID，值为几何对象
    source_geom_dict = nodes_gpd_ndh.set_index(source_id_field).geometry.to_dict()
    # 遍历目标shp文件，根据ID替换点坐标
    for idx, row in nodes_gpd.iterrows():
        target_id = row[target_id_field]
        if target_id in source_geom_dict:
            nodes_gpd.at[idx, 'geometry'] = source_geom_dict[target_id]
            nodes_gpd.at[idx, 'LGTD'] = source_geom_dict[target_id].x
            nodes_gpd.at[idx, 'LTTD'] = source_geom_dict[target_id].y
    # 覆盖原文件
    nodes_gpd.to_file(node_file, encoding="utf-8")
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
