"""
Author: Wenyu Ouyang
Date: 2025-02-01 08:18:13
LastEditTime: 2025-02-02 08:52:48
LastEditors: Wenyu Ouyang
Description: find the STCD of the nodes according to the index file
FilePath: \hydrotopo\scripts\find_node_id.py
Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
"""

import os
from pathlib import Path
import geopandas as gpd
from shapely.geometry import LineString
import pandas as pd

from hydroutils import hydro_file

# 设置文件路径
project_dir = Path(os.path.abspath(__file__)).parent.parent
data_dir = project_dir / "data"
result_dir = project_dir / "results"
sta_type = "RR"
node_file = (
    data_dir / f"{sta_type.lower()}_stations" / f"{sta_type.lower()}_stations.shp"
)
cur_sta = 2172
index_file = (
    result_dir / f"{cur_sta}_upstream_{sta_type}_station_lst.json"
)  # 假设索引文件名为 index_file.txt
output_file = result_dir / f"{cur_sta}_upstream_{sta_type}_stations_stcd.json"
line_shapefile = (
    project_dir / "results" / f"{cur_sta}_upstream_{sta_type}_station_connections"
)
point_shapefile = (
    project_dir / "results" / f"{cur_sta}_upstream_{sta_type}_station_points"
)

# 读取 GeoDataFrame
nodes_gpd = gpd.read_file(node_file)

# 读取索引文件
index_data = hydro_file.unserialize_json(index_file)

# 根据索引找到对应的 STCD
stcd_data_ = []
for index_list in index_data:
    stcd_list = [nodes_gpd.loc[idx, "STCD"] for idx in index_list]
    stcd_data_.append(stcd_list)
# find the name of the station according to the STCD-name mapping file
# 读取 STCD-Name 映射文件
stcd_name_file = data_dir / "ST_STBPRP_B.csv"
stcd_name_data = pd.read_csv(stcd_name_file, index_col=0)
stcd_name_dict = dict(zip(stcd_name_data["STCD"], stcd_name_data["STNM"]))
# 将 STCD 转换为站名
stcd_data = [
    [{stcd: stcd_name_dict[stcd].strip()} for stcd in stcd_list]
    for stcd_list in stcd_data_
]
# 将结果保存到文件
hydro_file.serialize_json(stcd_data, output_file, ensure_ascii=False)
print(f"STCD 数据已保存到 {output_file}")
# 每一行都是树的一个分支，我想根据其经纬度标出其位置，然后上下游之间连接成线
# 创建一个新的 GeoDataFrame 用于存储线条
lines_gdf = gpd.GeoDataFrame(columns=["geometry"], crs=nodes_gpd.crs)
points_gdf = gpd.GeoDataFrame(columns=["STCD", "NAME", "geometry"], crs=nodes_gpd.crs)

for branch in stcd_data:
    # 创建一个新的 GeoDataFrame 用于存储站点位置
    gdf = gpd.GeoDataFrame(columns=["STCD", "NAME", "geometry"])
    for stcd in branch:
        row = nodes_gpd[nodes_gpd["STCD"].str.strip() == list(stcd.keys())[0].strip()]
        if not row.empty:
            new_row = gpd.GeoDataFrame(
                [
                    {
                        "STCD": list(stcd.keys())[0].strip(),
                        "NAME": list(stcd.values())[0].strip(),
                        "geometry": row.geometry.values[0],
                    }
                ],
                crs=nodes_gpd.crs,
            )
            gdf = pd.concat([gdf, new_row], ignore_index=True)
            points_gdf = pd.concat([points_gdf, new_row], ignore_index=True)
    # 创建线条
    if len(gdf) > 1:
        line = LineString(gdf.geometry.tolist())
        new_line_gdf = gpd.GeoDataFrame([{"geometry": line}], crs=nodes_gpd.crs)
        lines_gdf = pd.concat([lines_gdf, new_line_gdf], ignore_index=True)

# delete duplicate points
points_gdf = points_gdf.drop_duplicates(subset=["STCD"])
if not line_shapefile.exists():
    # 保存线条到 shapefile
    lines_gdf.to_file(line_shapefile)
if not point_shapefile.exists():
    # 保存站点到 shapefile
    points_gdf.to_file(point_shapefile, encoding="utf-8")
print(f"线条已保存到 {line_shapefile}")
print(f"站点已保存到 {point_shapefile}")
