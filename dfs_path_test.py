import os
import sys

import bidict as bd
import matplotlib.pyplot as plt
import geopandas as gpd
import networkx as nx
import pandas as pd
from geopandas import GeoDataFrame
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points, split, snap
from shapely.wkt import loads

INPUT_NETWORK_FILE_SHP = os.path.relpath("test_data/near_changchun_cut.shp")  # network file, lines (no multiline objects)
# INPUT_VERTICES_FILE_SHP = os.path.relpath("near_changchun_vertices.shp")  # 从near_changchun_cut中抽取出的顶点
INPUT_NODE_FILE_SHP = os.path.relpath("test_data/near_changchun_dots.shp")


# END_NODE_TOLERANCE = 0.001  # NO TOLERANCE FOR LINES THAT ARE NOT SNAPPED!!!


def read_shape_file_to_graph(file_str):
    """Read shp file, make a geopandas dataframe and change the INFILE_LENGTH column name
     geometry options 'LineString', 'Point' """
    # make a geopandas dataframe from file
    gpd_df = gpd.read_file(file_str)
    nodeset = set()
    directed_graph = nx.DiGraph()  # 生成有向图
    for x in range(0, len(gpd_df)):
        u = gpd_df.geometry[x].coords[0]  # x起点和终点，-1代指最后一个索引
        nodeset.add(u)
        v = gpd_df.geometry[x].coords[-1]
        nodeset.add(v)
        if not (directed_graph.has_edge(u, v) | directed_graph.has_edge(v, u)):
            geom = gpd_df.geometry[x]
            directed_graph.add_edge(u, v)
            directed_graph.edges[u, v]['geometry'] = geom
    nodelist = list(nodeset)
    return directed_graph, nodelist


def geopandas_min_dist(point, gpd_dataframe, initial_buffer=500):
    """https://gis.stackexchange.com/questions/266730/filter-by-bounding-box-in-geopandas/266833
    Use the .cx index to find geoms inside the bounding box.
    Only use calculate the distance based on that.  Increase the bounding
    box size is you don't have a geometry.
    Return THE nearest geometry from a point"""
    buffer_steps = 500
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
    geoseries = gpd_df.loc[gpd_df['Dist'].idxmin()]  # 取最小距离所代表图元
    return geoseries


def gen_stations_graph(nodes_file_str, network_file_str):
    gpd_nodes_df = gpd.read_file(nodes_file_str)
    gpd_network_df = gpd.read_file(network_file_str)
    list_point_line_tuple = []
    for point in gpd_nodes_df.geometry:
        list_point_line_tuple.append((point.to_wkt(), geopandas_min_dist(point, gpd_network_df, 200).geometry.to_wkt()))
    graph_frame = gpd.GeoDataFrame(list_point_line_tuple, columns=['near_stations', 'nearest_line'])
    grouped_graph_frame: GeoDataFrame = gpd.GeoDataFrame(graph_frame.groupby(['nearest_line']))
    directed_graph = nx.DiGraph()
    for group in grouped_graph_frame:
        nearest_line: LineString = group['nearest_line'][0]
        for i in range(0, len(group['near_stations'])):
            new_line_string = snap(nearest_line, group['near_stations'][i],
                                   group['near_stations'][i].distance(nearest_line) + 0.0001)  # 待修改
            for line in gpd_network_df.geometry:
                if line == nearest_line:
                    line = new_line_string
    for x in range(0, len(gpd_network_df)):
        src = gpd_network_df.geometry[x].coords[0]  # x起点和终点，-1代指最后一个索引
        dest = gpd_network_df.geometry[x].coords[-1]
        stations_node_list = []
        for control_point in gpd_network_df.geometry[x].coords:
            if control_point in graph_frame['near_stations']:
                stations_node_list.append(control_point)
        directed_graph.add_edge(src, stations_node_list[0])
        for i in range(0, len(stations_node_list)):
            directed_graph.add_edge(stations_node_list[i], stations_node_list[i + 1])
        directed_graph.add_edge(stations_node_list[-1], dest)
    return directed_graph, grouped_graph_frame


'''def get_extrapoled_line(p1, p2, extrapol_ratio):  # 按照比率生成p1和p2之间的延长线
    """Creates a line extrapolated in p1->p2 direction，Input is shapely objects"""
    c = (p1.x + extrapol_ratio * (p2.x - p1.x), p1.y + extrapol_ratio * (p2.y - p1.y))
    return LineString([p1, c])'''

'''def tie_outside_node(gpd_df_nodes, gpd_df_network):
    """For each outside node, find the nearest line inside the buffer. The buffer
    grows if no lines are found.  """
    global nearest_point_exact
    vertices_df: GeoDataFrame = gpd.read_file(INPUT_VERTICES_FILE_SHP)  # 读入抽取出来的顶点文件
    for x in range(0, len(gpd_df_nodes)):
        source_target_nodes_geom = gpd_df_nodes.geometry[x]  # round geopandas rows'geometry to wkt
        nearest_line = geopandas_min_dist(gpd_df_nodes.geometry[x], gpd_df_network)
        nearest_line_geometry = nearest_line.geometry
        # 遍历顶点文件中的node,确定x,y是否都在其中
        if source_target_nodes_geom not in vertices_df.geometry:
            # Find the nearest point on this nearest line, find possible join_line.
            # This line will intersect only if there is a node on the line
            outside_node_grid_geom_to_nearest_p = \
                [q.wkt for q in nearest_points(source_target_nodes_geom, nearest_line_geometry)]
            # load point 2 (nearest point)
            nearest_p = loads(outside_node_grid_geom_to_nearest_p[1])
            # Make a coordinate tuple. Check if it's a first or last node
            near_first_node = nearest_p.distance((Point(nearest_line_geometry.coords[0])))
            near_last_node = nearest_p.distance((Point(nearest_line_geometry.coords[-1])))
            if (near_first_node < END_NODE_TOLERANCE) & (near_last_node >= END_NODE_TOLERANCE):
                nearest_point_exact = (Point(nearest_line_geometry.coords[0]))
            elif (near_first_node >= END_NODE_TOLERANCE) & (near_last_node < END_NODE_TOLERANCE):
                nearest_point_exact = (Point(nearest_line_geometry.coords[-1]))
            else:
                join_line = LineString([source_target_nodes_geom, nearest_p])
                splits_edges_coll = []
                start_buf_ratio = 1 + 0.001 / join_line.length  # grow 1 mm at a time
                # Break out of endless loop
                break_count = 0
                while (len(splits_edges_coll)) != 2:
                    join_line_interpolation = get_extrapoled_line(source_target_nodes_geom, nearest_p, start_buf_ratio)
                    splits_edges_coll = split(nearest_line_geometry, join_line_interpolation)
                    start_buf_ratio += 1 / join_line.length
                    break_count += 1
                    if break_count == 1000:
                        sys.exit('Did over 1000 loops.')
                nearest_point_exact = splits_edges_coll[0].coords[-1]  # Select last node in geom 0 of splits_edges_coll as
    return nearest_point_exact'''

'''def build_graph(file_str_network, file_str_node):  # 换成pd.dataframe
    old_nodes_df = gpd.read_file(file_str_node)
    old_edges_df = gpd.read_file(file_str_network)
    near_stations_from_src = []
    nearest_line_set = set()
    for i in range(0, len(old_nodes_df)):
        graph_frame['nearest_line'][i]: LineString = geopandas_min_dist(old_nodes_df[i], old_edges_df)
        graph_frame['source'][i] = graph_frame['nearest_line'][i].geometry[0]
        graph_frame['dest'][i] = graph_frame['nearest_line'][i].geometry[-1]
        # add near_stations
    for line in graph_frame['nearest_line']:
        nearest_line_set.add(line)  # 去重
    for line in nearest_line_set:
        nearest_line_set.add('')'''
'''
    network_graph: DiGraph = read_shape_file_to_graph(file_str_network)[0]
    source_nearest_dict = tie_outside_node(old_nodes_df, old_edges_df)[3]
    point_line_str_dict = tie_outside_node(old_nodes_df, old_edges_df)[2]
    list_of_line_geom = []  # 用来保存某条河段边上的所有点
    for source_point in source_nearest_dict.values():
        for nearest_point in source_nearest_dict.get(source_point):
            list_of_line_geom.append(nearest_point)
        list_of_line_geom.sort(key=lambda point: point.distance(source_point))
        for i in range(0, len(list_of_line_geom)):
            list_of_line_geom[i] = tie_outside_node(old_nodes_df, old_edges_df)[1].get(nearest_point)
        network_graph.add_edge(source_point, list_of_line_geom[0])
        for index in range(0, len(list_of_line_geom) - 2):
            network_graph.add_edge(list_of_line_geom[index], list_of_line_geom[index + 1])
        network_graph.add_edge(list_of_line_geom[-1], point_line_str_dict.get(nearest_point)[-1])
    return network_graph'''


def show_test_graph(graph):
    nx.draw(graph, node_size=10)
    plt.show()


def show_downstream_path(graph, node_list, index):
    down_path = nx.dfs_tree(graph, node_list[index])
    nx.draw(down_path, node_size=10)
    plt.show()
    path_coord_list = list(down_path)
    path_down_list = []
    for i in range(0, path_coord_list.__len__()):
        if path_coord_list[i] in node_list:
            path_down_list.append(node_list.index(path_coord_list[i]))
    print(path_down_list)


def show_upstream_paths(graph, node_list, index):
    path_up_list = []
    path_up_graph = graph.subgraph(nx.ancestors(graph, node_list[index]) | {node_list[index]}).copy()
    for coord in path_up_graph.nodes:
        if path_up_graph.in_degree(coord) == 0 & nx.generic.has_path(path_up_graph, coord, node_list[index]):
            for up_path in nx.all_simple_paths(path_up_graph, coord, node_list[index]):
                for path_coord in up_path:
                    path_up_list.append(node_list.index(path_coord))
                print(path_up_list)
                path_up_list.clear()
    nx.draw(path_up_graph, node_size=10)
    plt.show()
    print(path_up_list)


if __name__ == '__main__':
    # index = 27
    # test_graph = read_shape_file_to_graph(INPUT_NETWORK_FILE_SHP)[0]
    # test_node_list = read_shape_file_to_graph(INPUT_NETWORK_FILE_SHP)[1]
    # show_test_graph(test_graph)
    # show_downstream_path(test_graph, test_node_list, index)
    # show_upstream_paths(test_graph, test_node_list, index)
    test_graph1 = gen_stations_graph(INPUT_NODE_FILE_SHP, INPUT_NETWORK_FILE_SHP)
    show_test_graph(test_graph1)
    # test_graph.clear()
