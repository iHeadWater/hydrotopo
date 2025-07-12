"""
Author: Yang Wang
Date: 2024-12-12 17:03:14
LastEditTime: 2025-07-12 17:07:37
LastEditors: Wenyu Ouyang
Description: Create graph using hydrotopo
FilePath: /hydrotopo/hydrotopo/create_graph.py
Copyright (c) 2023-2026 Yang Wang. All rights reserved.
"""
import os
import networkx as nx
import geopandas as gpd
import numpy as np
from itertools import chain

import hydrotopo.ig_path as htip


def get_upstream_graph(basin_list, res_dir, network_path, node_path):
    """
    Create or load upstream graph for given basins.
    
    This function creates a directed graph representing upstream relationships
    between hydrological stations. It first checks if pre-computed graph files
    exist, and if not, generates the graph from network and node shapefiles.
    
    Parameters
    ----------
    basin_list : list
        List of basin station IDs for which to create the upstream graph.
    res_dir : str
        Directory path where graph files will be saved or loaded from.
    network_path : str
        Path to the network shapefile containing stream network data.
    node_path : str
        Path to the node shapefile containing station/node information.
        
    Returns
    -------
    tuple
        A tuple containing:
        - total_graph : networkx.DiGraph
            Directed graph representing upstream relationships between nodes.
        - basin_station_df : pandas.DataFrame
            DataFrame containing node information with columns:
            'node_id', 'station_id', 'basin_id', 'upstream_len'.
    """
    import pandas as pd
    nx_graph_path = os.path.join(res_dir, f'total_graph_{len(basin_list)}.gexf')
    basin_station_path = os.path.join(res_dir, f'basin_stations_{len(basin_list)}.csv')
    if (os.path.exists(nx_graph_path)) & (os.path.exists(basin_station_path)):
        total_graph = nx.read_gexf(nx_graph_path)
        basin_station_df = pd.read_csv(basin_station_path, engine='c')
        graph_tuple = (total_graph, basin_station_df)
    else:
        # 保存成csv文件，存储节点所对应的站名、流域、上游几个点
        total_graph = nx.DiGraph()
        index_station_dict = {}
        index_basin_dict = {}
        up_len_dict = {}
        network_features = gpd.read_file(network_path)
        node_features = gpd.read_file(node_path)
        basins = [sta_id.split('_')[-1] for sta_id in basin_list]
        upstream_graphs = prepare_graph(network_features, node_features, basins)
        for key in upstream_graphs.keys():
            upstream_graph = upstream_graphs[key]
            id_col = 'ID' if 'ID' in node_features.columns else 'STCD'
            # basin_id = node_features[id_col][node_features.index == key].values[0]
            if len(upstream_graph) == 0:
                upstream_graph = nodes_arr = [key]
            elif upstream_graph.dtype == 'O':
                    # upstream_graph = array(list1, list2, list3)
                nodes_arr = (
                    np.unique(list(chain.from_iterable(upstream_graph)))
                    if upstream_graph.shape[0] > 1
                    else upstream_graph
                )
            else:
                # upstream_graph = array(list1, list2) and dtype is not object
                nodes_arr = np.unique(upstream_graph)
            nodes_arr = np.append(nodes_arr, key) if key not in nodes_arr else nodes_arr
            for node in nodes_arr:
                node_id = node_features[id_col][node_features.index == node].values[0]
                basin = node_id if node_id in basins else node_features[id_col][node_features.index == key].values[0]
                index_station_dict[int(node)] = (
                    f'songliao_{node_id}' if '_' not in node_id and f'songliao_{node_id}' in basin_list
                    else f'camels_{node_id}' if '_' not in node_id
                    else node_id)
                index_basin_dict[int(node)] = f'songliao_{basin}' if f'songliao_{basin}' in basin_list else f'camels_{basin}'
            up_len_dict[int(key)] = len(nodes_arr)
            for path in upstream_graph:
                path = [path] if isinstance(path, int | np.int64) else path
                path = np.append(path, key) if key not in path else path
                nx.add_path(total_graph, path)
        basin_station_df = pd.DataFrame([index_station_dict, index_basin_dict, up_len_dict]).T
        basin_station_df = basin_station_df[~basin_station_df[0].isna()].rename(columns={0: 'station_id', 1: 'basin_id', 2: 'upstream_len'})
        basin_station_df = basin_station_df.reset_index().rename(columns={'index': 'node_id'})
        graph_tuple = (total_graph, basin_station_df)
        # 有孤立点的情况不适用于edgelist, int点不适合gml
        nx.write_gexf(total_graph, nx_graph_path)
        basin_station_df.to_csv(basin_station_path)
    return graph_tuple

def prepare_graph(network_features: gpd.GeoDataFrame, node_features: gpd.GeoDataFrame, nodes: list[int] | list[str],
                  cutoff=4):
    """
    Prepare graph data structure for upstream path analysis.
    
    This function processes network and node geographic data to create a graph
    dictionary containing upstream relationships. It handles both integer and
    string node identifiers and uses the hydrotopo.ig_path module for bulk
    upstream path finding.
    
    Parameters
    ----------
    network_features : geopandas.GeoDataFrame
        GeoDataFrame containing network/stream line features with geometric
        information for building the graph structure.
    node_features : geopandas.GeoDataFrame
        GeoDataFrame containing node/station point features with ID information.
        Must contain either 'ID' or 'STCD' column for node identification.
    nodes : list of int or list of str
        List of node identifiers to process. Can be either integer indices
        or string IDs that will be matched against the node_features DataFrame.
    cutoff : int, optional
        Maximum depth for upstream path searching. Default is 4.
        
    Returns
    -------
    dict
        Dictionary mapping node indices to their upstream path information.
        Keys are node indices, values contain upstream connectivity data
        as returned by hydrotopo.ig_path.find_edge_nodes_bulk_up.
    """
    # test_df_path = 's3://stations-origin/zq_stations/hour_data/1h/zq_CHN_songliao_10800300.csv'
    if isinstance(nodes[0], int):
        node_idx = nodes
    else:
        id_col = 'ID' if 'ID' in node_features.columns else 'STCD'
        node_features[id_col] = node_features[id_col].astype(str)
        node_idx = node_features.index[node_features[id_col].isin(nodes)]
    return (
        htip.find_edge_nodes_bulk_up(
            node_features, network_features, node_idx, cutoff
        )
        if len(node_idx) != 0
        else {}
    )
