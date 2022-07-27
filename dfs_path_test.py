import os
import sys

import matplotlib.pyplot as plt
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points, split
from shapely.wkt import loads

INPUT_NETWORK_FILE_SHP = os.path.relpath("Sample_Network3D.shp")  # network file, lines (no multiline objects)
INPUT_NODE_FILE_SHP = os.path.relpath("near_changchun_dots.shp")
END_NODE_TOLERANCE = 0.01  # Tolerance for end node snap around line 490 only NO TOLERANCE FOR LINES THAT ARE NOT SNAPPED!!!


def read_shape_file_to_graph(file_str):
    """Read shp file, make a geopandas dataframe
    and change the INFILE_LENGTH column name
    geometry options 'LineString', 'Point'
    """
    # make a geopandas dataframe from file
    gpd_df = gpd.read_file(file_str)
    nodeset = set()
    directed_graph = nx.DiGraph()  # 生成有向图
    for x in range(0, len(gpd_df)):
        u = gpd_df.geometry[x].coords[0]  # x起点和终点，-1代指最后一个索引
        nodeset.add(u)
        v = gpd_df.geometry[x].coords[-1]
        nodeset.add(v)
        if directed_graph.has_edge(u, v) | directed_graph.has_edge(v, u):
            pass
        else:
            geom = gpd_df.geometry[x]
            directed_graph.add_edge(u, v)
            directed_graph.edges[u, v]['geometry'] = geom
    nodelist = list(nodeset)
    return directed_graph, nodelist


def geopandas_min_dist(point, gpd_dataframe, initial_buffer=1000):
    """https://gis.stackexchange.com/questions/266730/filter-by-bounding-box-in-geopandas/266833
    Use the .cx index to find geoms inside the bounding box.
    Only use calculate the distance based on that.  Increase the bounding
    box size is you don't have a geometry.
    Return THE nearest geometry from point a point
    """
    buffer_steps = 1000
    # 给点或者其他geometry加buffer
    xmin, ymin, xmax, ymax = (point.buffer(initial_buffer)).bounds
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
    geoseries = gpd_df.loc[gpd_df['Dist'].idxmin()]  # 取最小距离
    return geoseries


def get_extrapoled_line(p1, p2, extrapol_ratio):  # 按照比率生成p1和p2之间的截线p1c
    """Creates a line extrapolated in p1->p2 direction，Input is shapely objects"""
    c = (p1.x + extrapol_ratio * (p2.x - p1.x), p1.y + extrapol_ratio * (p2.y - p1.y))
    return LineString([p1, c])


def tie_outside_node(gpd_df_nodes, gpd_df_network):
    """For each outside node, find the nearest line inside the buffer. The buffer
    grows if no lines are found.  """
    graph = read_shape_file_to_graph(INPUT_NETWORK_FILE_SHP)
    for x in range(0, len(gpd_df_nodes)):
        source_target_nodes_geom = gpd_df_nodes.geometry[x]  # round geopandas rows'geometry to wkt
        nearest_line = geopandas_min_dist(gpd_df_nodes.geometry[x], gpd_df_network, 1000)
        nearest_line_geometry = nearest_line.geometry
        # 遍历graph中的node,确定x,y是否都在其中
        if (source_target_nodes_geom.x, source_target_nodes_geom.y) not in list(graph):
            # Find the nearest point on this nearest line, find possible join_line.
            # This line will intersect only if there is a node on the line
            outside_node_grid_geom_to_nearest_p = \
                [q.wkt for q in nearest_points(source_target_nodes_geom, nearest_line_geometry)]
            # load point 2 (nearest point)
            nearest_p = loads(outside_node_grid_geom_to_nearest_p[1])
            # Make a coordinate tuple. Check if it's a first or last node
            near_first_node = nearest_p.distance((Point(nearest_line_geometry.coords[0])))
            near_last_node = nearest_p.distance((Point(nearest_line_geometry.coords[-1])))
            if near_first_node <= near_last_node:
                nearest_point_case1 = (Point(nearest_line_geometry.coords[0]))
            else:
                nearest_point_case1 = (Point(nearest_line_geometry.coords[-1]))
            graph.add_edge(source_target_nodes_geom, nearest_point_case1)
            if not (near_first_node <= END_NODE_TOLERANCE or near_last_node <= END_NODE_TOLERANCE):
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
                nearest_point_exact_case2 = splits_edges_coll[0].coords[-1]  # Select last node in geom 0 of splits_edges_coll as
                graph.add_edge(source_target_nodes_geom, nearest_point_exact_case2)
    return graph


if __name__ == '__main__':
    # gpddf_nodes = gpd.read_file(INPUT_NODE_FILE_SHP)
    # gpddf_network = gpd.read_file(INPUT_NETWORK_FILE_SHP)
    # tie_outside_node(gpddf_nodes, gpddf_network)
    test_graph = read_shape_file_to_graph(INPUT_NETWORK_FILE_SHP)[0]
    nx.draw(test_graph, node_size=5)
    plt.show()
    index = 27
    test_node_list = read_shape_file_to_graph(INPUT_NETWORK_FILE_SHP)[1]
    path = nx.dfs_tree(test_graph, test_node_list[index])
    nx.draw(path, node_size=5)
    plt.show()
    path_coord_list = list(path)
    path_down_list = []
    for i in range(0, path_coord_list.__len__()):
        if path_coord_list[i] in test_node_list:
            path_down_list.append(test_node_list.index(path_coord_list[i]))
    print(path_down_list)
    for i in range(0, test_graph.number_of_nodes()):
        if nx.generic.has_path(test_graph, test_node_list[i], test_node_list[index]) & i != index:
            path_up_list = []
            for path in nx.all_simple_paths(test_graph, test_node_list[i], test_node_list[index]):
                for coord in path:
                    path_up_list.append(test_node_list.index(coord))
                print(path_up_list)
                path_up_list.clear()
    path_up_list = []
    path_up_graph = test_graph.subgraph(nx.ancestors(test_graph, test_node_list[index])).copy()
    for coord in path_up_graph.nodes:
        if path_up_graph.in_degree(coord) == 0 & nx.generic.has_path(path_up_graph, coord, test_node_list[index]):
            for up_path in nx.all_simple_paths(path_up_graph, coord, test_node_list[index]):
                for path_coord in up_path:
                    path_up_list.append(test_node_list.index(path_coord))
                    print(path_up_list)
                    path_up_list.clear()
    nx.draw(path_up_graph, node_size=5)
    plt.show()
    print(path_up_list)
    test_graph.clear()
