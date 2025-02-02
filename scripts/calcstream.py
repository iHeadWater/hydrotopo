"""
Author: Wenyu Ouyang
Date: 2025-01-28 12:27:59
LastEditTime: 2025-02-03 06:29:40
LastEditors: Wenyu Ouyang
Description: try to use cli.py to run the function
FilePath: \hydrotopo\scripts\calcstream.py
Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
"""

import os
from pathlib import Path
import sys
import geopandas as gpd
import pandas as pd

from hydroutils.hydro_file import serialize_json_np

sys.path.append(os.path.dirname(Path(os.path.abspath(__file__)).parent))
from hydrotopo.ig_path import find_main_and_tributary, find_edge_nodes


def calcstream(
    nodes_path, river_path, cur_sta, up_sta, cutoff, sta_type, upstream, downstream
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
        result_dir = project_dir / "results"
        save_json_file = result_dir / f"{cur_sta}_upstream_{sta_type}_station_lst.json"
        if not save_json_file.exists():
            serialize_json_np(upstream_station_lst, save_json_file)
        print(upstream_station_lst)
    elif downstream is True:
        print(find_edge_nodes(nodes_gpd, network_gpd, cur_sta, "down", cutoff))
    if up_sta is not None:
        print(find_main_and_tributary(nodes_gpd, network_gpd, cur_sta, up_sta))


if __name__ == "__main__":
    project_dir = Path(os.path.abspath(__file__)).parent.parent
    data_dir = project_dir / "data"
    result_dir = project_dir / "results"
    cur_sta_lst = [2172]
    sta_type_lst = ["RR"]
    river_file = result_dir / "northeast_rivers" / "northeast_rivers.shp"
    for cur_sta in cur_sta_lst:
        for sta_type in sta_type_lst:
            node_file = (
                data_dir
                / f"{sta_type.lower()}_stations"
                / f"{sta_type.lower()}_stations.shp"
            )
            calcstream(
                nodes_path=node_file,
                river_path=river_file,
                cur_sta=cur_sta,
                up_sta=None,
                cutoff=2,
                sta_type=sta_type,
                upstream=True,
                downstream=False,
            )
