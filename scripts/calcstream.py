"""
Author: Wenyu Ouyang
Date: 2025-01-28 12:27:59
LastEditTime: 2025-02-04 15:35:29
LastEditors: Wenyu Ouyang
Description: try to use cli.py to run the function
FilePath: \hydrotopo\scripts\calcstream.py
Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
"""

import os
import itertools
import sys
import pandas as pd
import geopandas as gpd
from pathlib import Path

from hydroutils.hydro_file import serialize_json_np

sys.path.append(os.path.dirname(Path(os.path.abspath(__file__)).parent))
from hydrotopo.ig_path import find_main_and_tributary, find_edge_nodes


def calcstream(
    nodes_path,
    river_path,
    cur_sta,
    up_sta,
    cutoff,
    upstream,
    downstream,
    up_save_file=None,
    down_save_file=None,
):
    """calculate the upstream or downstream stations of the given station

    Parameters
    ----------
    nodes_path : str
        the path of the nodes shape file
    river_path : str
        the path of the river shape
    cur_sta : int
        the index of current station
    up_sta : int
        the index of the upstream station
        if given, find the path between the current station and the upstream station
    cutoff : int
        the maximum distance (number of nodes) between the current station and the target station
        including the current station
    sta_type : str
        RR, ZZ, ZQ or others, just for naming the output file now
    upstream : bool
        if True, calculate the upstream stations
    downstream : bool
        if True, calculate the downstream stations
    up_save_file : str
        the path to save the upstream stations
    down_save_file : str
        the path to save the downstream stations

    Raises
    ------
    ValueError
        the input shape file contains nan values
    ValueError
        the input shape file contains nan values
    """
    input_network_file_shp = os.path.abspath(river_path)
    input_node_file_shp = os.path.relpath(nodes_path)
    nodes_gpd = gpd.read_file(input_node_file_shp)
    network_gpd = gpd.read_file(input_network_file_shp)
    # check if any line is nan
    if nodes_gpd.isnull().values.any():
        raise ValueError("There are nan values in nodes shape file; please delete them")
    if network_gpd.isnull().values.any():
        raise ValueError("There are nan values in river shape file; please delete them")
    if upstream is True:
        upstream_station_lst = find_edge_nodes(
            nodes_gpd, network_gpd, cur_sta, "up", cutoff
        )
        # save the result to a file
        if up_save_file is not None and not up_save_file.exists():
            serialize_json_np(upstream_station_lst, up_save_file)
        print(upstream_station_lst)
    elif downstream is True:
        print(find_edge_nodes(nodes_gpd, network_gpd, cur_sta, "down", cutoff))
    if up_sta is not None:
        print(find_main_and_tributary(nodes_gpd, network_gpd, cur_sta, up_sta))


if __name__ == "__main__":
    project_dir = Path(os.path.abspath(__file__)).parent.parent
    data_dir = project_dir / "data"
    result_dir = project_dir / "results"
    sta_type_lst = ["RR", "ZZ"]
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
    river_file = result_dir / "northeast_rivers" / "northeast_rivers.shp"
    new_node_dir = result_dir / "tmp_nodes"
    new_node_dir.mkdir(parents=True, exist_ok=True)
    for stcd, sta_type in itertools.product(target_stcd_lst, sta_type_lst):
        node_file = (
            data_dir
            / f"{sta_type.lower()}_stations"
            / f"{sta_type.lower()}_stations.shp"
        )
        nodes_gpd = gpd.read_file(node_file)
        # find the row index of the target node with its STCD
        cur_sta = nodes_gpd[nodes_gpd["STCD"] == stcd].index
        new_node_file = new_node_dir / f"{stcd}_{sta_type}_stations.shp"
        # if the station is not found, get its info from target_stations.shp and add it to the stations shp
        if cur_sta is None or len(cur_sta) == 0:
            target_nodes_gpd = gpd.read_file(result_dir / "target_stations")
            target_node = target_nodes_gpd[target_nodes_gpd["STCD"] == stcd]
            nodes_gpd_new = pd.concat([nodes_gpd, target_node], ignore_index=True)
            nodes_gpd_new.to_file(new_node_file, encoding="utf-8")
            # update cur_sta, use the index of the new station
            cur_sta = nodes_gpd_new[nodes_gpd_new["STCD"] == stcd].index[0]
        else:
            # copy the file to the new one
            nodes_gpd.to_file(new_node_file, encoding="utf-8")
            cur_sta = cur_sta[0]
        up_save_file = result_dir / f"{stcd}_upstream_{sta_type}_station_lst.json"
        calcstream(
            nodes_path=new_node_file,
            river_path=river_file,
            cur_sta=cur_sta,
            up_sta=None,
            cutoff=2,
            upstream=True,
            downstream=False,
            up_save_file=up_save_file,
        )
