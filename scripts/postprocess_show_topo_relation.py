"""
Author: Wenyu Ouyang
Date: 2025-02-01 08:18:13
LastEditTime: 2025-02-04 18:30:15
LastEditors: Wenyu Ouyang
Description: find the STCD of the nodes according to the index file
FilePath: \hydrotopo\scripts\postprocess_show_topo_relation.py
Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
"""

import os
import itertools
import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import LineString

from hydroutils import hydro_file

# 设置文件路径
project_dir = Path(os.path.abspath(__file__)).parent.parent
data_dir = project_dir / "data"
result_dir = project_dir / "results"


def postprocess(stcd="20600340", sta_type="RR"):
    tmp_node_dir = result_dir / "tmp_nodes"
    node_file = tmp_node_dir / f"{stcd}_{sta_type}_stations.shp"
    index_file = (
        result_dir / f"{stcd}_upstream_{sta_type}_station_lst.json"
    )  # 假设索引文件名为 index_file.txt
    output_file = result_dir / f"{stcd}_upstream_{sta_type}_stations_stcd.json"
    line_shapefile = (
        project_dir / "results" / f"{stcd}_upstream_{sta_type}_station_connections"
    )
    point_shapefile = (
        project_dir / "results" / f"{stcd}_upstream_{sta_type}_station_points"
    )

    # 读取 GeoDataFrame
    nodes_gpd = gpd.read_file(node_file)

    # 读取索引文件
    index_data = hydro_file.unserialize_json(index_file)

    # 根据索引找到对应的 STCD
    stcd_data_ = []
    for index_list in index_data:
        if not isinstance(index_list, list):
            # sometimes the index_list is not a list only a single number when there is only one station
            index_list = [index_list]
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
    points_gdf = gpd.GeoDataFrame(
        columns=["STCD", "NAME", "geometry"], crs=nodes_gpd.crs
    )

    for branch in stcd_data:
        # 创建一个新的 GeoDataFrame 用于存储站点位置
        gdf = gpd.GeoDataFrame(columns=["STCD", "NAME", "geometry"])
        for stcd in branch:
            row = nodes_gpd[
                nodes_gpd["STCD"].str.strip() == list(stcd.keys())[0].strip()
            ]
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
    # 保存线条到 shapefile
    lines_gdf.to_file(line_shapefile)
    # 保存站点到 shapefile
    points_gdf.to_file(point_shapefile, encoding="utf-8")
    print(f"线条已保存到 {line_shapefile}")
    print(f"站点已保存到 {point_shapefile}")


if __name__ == "__main__":
    sta_type_lst = ["RR", "ZZ", "ZQ"]
    # 9 main reservoirs in Liaoning: 石佛寺，柴河，清河，闹德海，大伙房，观音阁，葠窝水库，汤河水库，白石水库
    target_stcd_lst = [
        # "20600340",
        # "20800900",
        # "20810200",
        "20910930",
        # "21100150",
        # "21110150",
        # "21110400",
        # "21113800",
        # "21200510",
    ]
    for stcd, sta_type in itertools.product(target_stcd_lst, sta_type_lst):
        postprocess(stcd=stcd, sta_type=sta_type)
    print("postprocess_find_node_id.py has been executed")
