import os
import sys

import bidict as bd
import geopandas as gpd
import matplotlib.pyplot as plt
import multidict as md
import networkx as nx
from networkx import DiGraph
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points, split
from shapely.wkt import loads

# 本文件测试数据来自https://www.hydrosheds.org/products/hydrorivers，坐标为EPSG:4326，单位全部为度，1°≈110km
INPUT_NETWORK_FILE_SHP = os.path.relpath("test_data/near_changchun_cut.shp")
INPUT_NODE_FILE_SHP = os.path.relpath("test_data/near_changchun_dots.shp")


def geopandas_min_dist(point, gpd_dataframe, initial_buffer=0.005):
    """https://gis.stackexchange.com/questions/266730/filter-by-bounding-box-in-geopandas/266833
    从限定范围内获取图元，本例中为河网线"""
    buffer_steps = 0.002
    # 给点或者其他geometry加buffer
    xmin, ymin, xmax, ymax = point.buffer(initial_buffer).bounds
    # 使用gpd_dataframe.cx获得特定buffer中的所有geometry
    gpd_df = gpd_dataframe.cx[xmin:xmax, ymin:ymax]
    # if empty, go into loop
    while gpd_df.empty:
        initial_buffer = initial_buffer + buffer_steps
        xmin, ymin, xmax, ymax = (point.buffer(initial_buffer)).bounds
        gpd_df = gpd_dataframe.cx[xmin:xmax, ymin:ymax]
    gpd_df = gpd_df.copy(deep=True)
    # 将每行中的geometry提取出来，与待测点做距离，添加到新列[Dist]中，axis=1代表对列做出修改
    gpd_df['Dist'] = gpd_df.apply(lambda row: point.distance(row.geometry), axis=1)
    geoseries = gpd_df.loc[gpd_df['Dist'].idxmin()]  # 取最小距离所代表图元，为dataframe中一整行
    return geoseries


def get_extrapolated_line(p1, p2, extrapolate_ratio):  # 按照比率生成p1和p2之间的延长线
    c = (p1.x + extrapolate_ratio * (p2.x - p1.x), p1.y + extrapolate_ratio * (p2.y - p1.y))
    return LineString([p1, c])


def tie_outside_node(gpd_df_nodes, gpd_df_network):  # 将外部站点绑定到河网线上
    source_point_line_dict = bd.bidict()  # 构建原点和投影点之间的映射关系
    nearest_point_line_dict = md.MultiDict()  # 构建投影点和线的映射关系
    for x in range(0, len(gpd_df_nodes)):
        source_target_nodes_geom = gpd_df_nodes.geometry[x]
        nearest_line = geopandas_min_dist(gpd_df_nodes.geometry[x], gpd_df_network).geometry
        # 用以解决低质数据，下文while处若循环1000次后还连不上，就放弃
        outside_node_grid_geom_to_nearest_p \
            = [q.wkt for q in nearest_points(source_target_nodes_geom, nearest_line)]
        # 从两个最近点中取一个
        nearest_p = loads(outside_node_grid_geom_to_nearest_p[1])
        join_line = LineString([source_target_nodes_geom, nearest_p])
        splits_edges_coll = []
        start_buf_ratio = 1 + 0.00001 / join_line.length
        break_count = 0
        while (len(splits_edges_coll)) != 2:
            join_line_interpolation = get_extrapolated_line(source_target_nodes_geom, nearest_p, start_buf_ratio)
            splits_edges_coll = split(nearest_line, join_line_interpolation)
            start_buf_ratio += 1 / join_line.length
            break_count += 1
            if break_count == 1000:
                sys.exit('Did over 1000 loops.')
        nearest_point_exact = splits_edges_coll[0].coords[-1]
        source_point_line_dict.put((source_target_nodes_geom.x, source_target_nodes_geom.y), nearest_point_exact)
        nearest_point_line_dict.add(nearest_line.to_wkt(), nearest_point_exact)
    return nearest_point_line_dict, source_point_line_dict


def build_graph(file_str_network, file_str_node):
    nodes_df = gpd.read_file(file_str_node)
    edges_df = gpd.read_file(file_str_network)
    network_graph: DiGraph = nx.DiGraph()
    nearest_point_line_dict = tie_outside_node(nodes_df, edges_df)[0]
    source_point_line_dict = tie_outside_node(nodes_df, edges_df)[1]
    for line in edges_df.geometry:
        line_wkt = line.to_wkt()
        list_point_and_dis = []
        list_coords = []
        if line_wkt in nearest_point_line_dict.keys():
            for coord in nearest_point_line_dict.getall(line_wkt):
                list_point_and_dis.append((coord, Point(coord).distance(Point(line.coords[0]))))
            list_point_and_dis.append((line.coords[0], 0))
            list_point_and_dis.append((line.coords[-1], 360))
            list_point_and_dis = sorted(list_point_and_dis, key=lambda point_and_dis: point_and_dis[1])
            for i in range(0, len(list_point_and_dis)):
                list_coords.append(list_point_and_dis[i][0])
            nx.add_path(network_graph, list_coords)
        else:
            src = line.coords[0]
            dest = line.coords[-1]
            network_graph.add_edge(src, dest)
    return network_graph, source_point_line_dict


def get_upstream_stations(nodes_file_path, network_file_path, station_index: int, cutoff: int = 2147483547):
    stations_up_list = []
    stations_graph = build_graph(network_file_path, nodes_file_path)[0]
    gpd_nodes_df = gpd.read_file(nodes_file_path)
    source_target_dict = build_graph(network_file_path, nodes_file_path)[1]
    source_node_coord = (gpd_nodes_df.geometry[station_index].x, gpd_nodes_df.geometry[station_index].y)
    target_node_coord = source_target_dict.get(source_node_coord)
    upstream_graph = stations_graph.subgraph(nx.ancestors(stations_graph, target_node_coord) | {target_node_coord}).copy()
    set_up_no_dup = set()  # 给找到的路径去重
    for coord in upstream_graph.nodes:
        if upstream_graph.in_degree(coord) == 0 & nx.generic.has_path(upstream_graph, coord, target_node_coord):
            for up_path in nx.all_simple_paths(upstream_graph, coord, target_node_coord):
                for path_coord in up_path:
                    if path_coord in source_target_dict.values():
                        point = Point(source_target_dict.inv[path_coord])
                        for x in range(0, len(gpd_nodes_df)):
                            if gpd_nodes_df.geometry[x] == point:
                                stations_up_list.append('station' + str(x + 1))
                if len(stations_up_list) > 1:
                    set_up_no_dup.add(str(stations_up_list[-cutoff:]))
                stations_up_list.clear()
    return upstream_graph, set_up_no_dup


def show_downstream_stations(nodes_file_path, network_file_path, station_index: int, cutoff: int = 2147483647):
    stations_graph = build_graph(network_file_path, nodes_file_path)[0]
    source_target_dict = build_graph(network_file_path, nodes_file_path)[1]
    gpd_nodes_df = gpd.read_file(nodes_file_path)
    source_node_coord = (gpd_nodes_df.geometry[station_index].x, gpd_nodes_df.geometry[station_index].y)
    target_node_coord = source_target_dict.get(source_node_coord)
    dfs_tree_path = nx.dfs_tree(stations_graph, target_node_coord)
    list_stations = []
    for coord in dfs_tree_path:
        if coord in source_target_dict.values():
            source_coord = source_target_dict.inv[coord]
            point = Point(source_coord)
            for x in range(0, len(gpd_nodes_df)):
                if gpd_nodes_df.geometry[x] == point:
                    list_stations.append('station' + str(x + 1))
                    break
    print(list_stations[:cutoff])


def show_upstream_stations_graph(nodes_file_path, network_file_path, number: int, cutoff: int):
    upstream_graph = get_upstream_stations(nodes_file_path, network_file_path, number, cutoff)[0]
    set_up_no_dup = get_upstream_stations(nodes_file_path, network_file_path, number, cutoff)[1]
    origin_graph = nx.DiGraph()
    new_graph = nx.DiGraph()
    gpd_df_lines = gpd.read_file(network_file_path)
    gpd_nodes_df = gpd.read_file(nodes_file_path)
    for line in gpd_df_lines.geometry:
        src = line.coords[0]
        dest = line.coords[-1]
        origin_graph.add_edge(src, dest)
    source_target_dict = build_graph(network_file_path, nodes_file_path)[1]
    source_node_coord = (gpd_nodes_df.geometry[number].x, gpd_nodes_df.geometry[number].y)
    target_node_coord = source_target_dict.get(source_node_coord)
    for coord in upstream_graph.nodes:
        if upstream_graph.in_degree(coord) == 0:
            for up_path in nx.all_simple_paths(upstream_graph, coord, target_node_coord):
                list_up_path = list(up_path)
                for path_coord in up_path:
                    if path_coord in origin_graph.nodes:
                        list_up_path.remove(path_coord)
                nx.add_path(new_graph, list_up_path[-cutoff:])
    for list_str in set_up_no_dup:
        print(list_str.lstrip('[').rstrip(']'))   # 输出的都是字串
    nx.draw(new_graph, node_size=10)
    plt.show()


if __name__ == '__main__':
    index = 0
    show_upstream_stations_graph(INPUT_NODE_FILE_SHP, INPUT_NETWORK_FILE_SHP, index, 3)
    show_downstream_stations(INPUT_NODE_FILE_SHP, INPUT_NETWORK_FILE_SHP, index, 4)

