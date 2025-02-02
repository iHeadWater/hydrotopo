"""
Author: Wenyu Ouyang
Date: 2025-01-28 12:27:59
LastEditTime: 2025-02-02 08:44:28
LastEditors: Wenyu Ouyang
Description: try to use cli.py to run the function
FilePath: \hydrotopo\scripts\calcstream.py
Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
"""

import os
from pathlib import Path
import sys
import click
import geopandas as gpd

from hydroutils.hydro_file import serialize_json_np

sys.path.append(os.path.dirname(Path(os.path.abspath(__file__)).parent))
from hydrotopo.ig_path import find_main_and_tributary, find_edge_nodes

project_dir = Path(os.path.abspath(__file__)).parent.parent
data_dir = project_dir / "data"
result_dir = project_dir / "results"
cur_sta_index = 2172
sta_type = "RR"
node_file = (
    data_dir / f"{sta_type.lower()}_stations" / f"{sta_type.lower()}_stations.shp"
)
river_file = result_dir / "northeast_rivers" / "northeast_rivers.shp"


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--nodes_path", help="path of nodes shape file", default=node_file)
@click.option(
    "--river_path", help="path of river vector shape file", default=river_file
)
@click.option(
    "--cur_sta", type=int, help="Index of current station", default=cur_sta_index
)
@click.option(
    "--up_sta",
    type=int,
    help="number of station which will be judge in mainstream or tributary in upstream "
    "watershed of current station",
    default=None,
)
@click.option(
    "--cutoff",
    type=int,
    default=2,
    help="amount of stations which user want to limit, including cur_sta itself",
)
@click.option(
    "--upstream",
    default=True,
    help="output upstream stations graph of current station",
)
@click.option(
    "--downstream",
    default=False,
    help="output list of downstream stations of current station",
)
def main(nodes_path, river_path, cur_sta, up_sta, cutoff, upstream, downstream):
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
        click.echo(str(upstream_station_lst))
    elif downstream is True:
        click.echo(
            str(find_edge_nodes(nodes_gpd, network_gpd, cur_sta, "down", cutoff))
        )
    if up_sta is not None:
        click.echo(
            str(find_main_and_tributary(nodes_gpd, network_gpd, cur_sta, up_sta))
        )


if __name__ == "__main__":
    main()
