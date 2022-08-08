import os
import sys
import time

import bidict as bd
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import multidict as md
import networkx as nx
from geopandas import GeoDataFrame
from networkx import DiGraph
from shapely import geometry
from shapely.geometry import Point, LineString, GeometryCollection
from shapely.ops import nearest_points, split
from shapely.wkt import loads

"""
本文件线元测试数据来自https://www.hydrosheds.org/products/hydrorivers，坐标为EPSG:4326，单位全部为度，1°≈110km
河网图层只接受LineString，不接受MultiLineString，输入前请将其转化为单部件
"""

INPUT_NETWORK_FILE_SHP = os.path.relpath("test_data/near_changchun_cut.shp")
INPUT_NODE_FILE_SHP = os.path.relpath("test_data/near_changchun_dots.shp")
gpd_nodes_dataframe = gpd.read_file(INPUT_NODE_FILE_SHP)
gpd_network_dataframe = gpd.read_file(INPUT_NETWORK_FILE_SHP)
'''用以标记数据是否过期，默认为False，改为True并删除所有缓存文件后，将从头生成所有数据'''
outdated = False


def geopandas_min_dist(point, gpd_dataframe, initial_buffer=0.005):
    """
    https://gis.stackexchange.com/questions/266730/filter-by-bounding-box-in-geopandas/266833
    从限定范围内获取图元，本例中为河网线

    参数
    ----------
    point : shapely.POINT
        待寻找最近线段的点
    gpd_dataframe : geopandas.GeoDataFrame
        存储河网线段的GeoDataFrame数据结构
    initial_buffer : float
        给点加的buffer（搜寻半径）大小，初始为0.005度，约为550米

    返回
    ----------
    GeoDataFrame中的一个Series（对应图层属性表中一行+geometry项），实际运行中通常取其geometry项
    """
    buffer_steps = 0.002
    # 给点或者其他geometry加buffer
    xmin, ymin, xmax, ymax = point.buffer(initial_buffer).bounds
    # 使用gpd_dataframe.cx获得特定buffer中的所有geometry
    gpd_df = gpd_dataframe.cx[xmin:xmax, ymin:ymax]
    while gpd_df.empty:
        initial_buffer = initial_buffer + buffer_steps
        xmin, ymin, xmax, ymax = (point.buffer(initial_buffer)).bounds
        gpd_df = gpd_dataframe.cx[xmin:xmax, ymin:ymax]
    gpd_df = gpd_df.copy(deep=True)
    # 将每行中的geometry提取出来，与待测点做距离，添加到新列[Dist]中，axis=1代表对列做出修改
    gpd_df['Dist'] = gpd_df.apply(lambda row: point.distance(row.geometry), axis=1)
    # 取最小距离所代表图元，为dataframe中一整行
    geoseries = gpd_df.loc[gpd_df['Dist'].idxmin()]
    return geoseries


def get_extrapolated_line(source_coord, coord, extrapolate_ratio):
    """
    参数
    ----------
    source_coord: 源点坐标，为浮点2元组（tuple）
        待寻找最近线段的点
    coord: 另一点坐标，也是浮点2元组
    extrapolate_ratio: 延长率

    返回
    ----------
    由source_coord和coord确定的一条线段延长线
    """
    source = Point(source_coord)
    p2 = Point(coord)
    c = (source.x + extrapolate_ratio * (p2.x - source.x), source.y + extrapolate_ratio * (p2.y - source.y))
    return LineString([source, c])


def tie_outside_node(gpd_df_nodes, gpd_df_network):
    """
    值得注意，由于tie_outside_node操作很费时，故工程中如有source_project_points.csv文件和nearest_line_project_points.csv文件，
    程序会优先从这两个文件中读取数据，所以如要从头生成，必须删除这两个文件
    参数
    ----------
    gpd_df_nodes: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    gpd_df_network: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe

    返回
    ----------
    nearest_point_line_dict: 线上投影点和最近河流线的对应列表
    source_point_line_dict: 原站点和站点在线上投影点的对应列表
    """
    if (outdated is False) | ((outdated is True) & os.path.exists('source_project_points.csv') & os.path.exists(
            'nearest_line_project_points.csv')):
        source_point_frame = pd.read_csv('source_project_points.csv')
        line_point_frame = pd.read_csv('nearest_line_project_points.csv')
        source_point_line_dict = bd.bidict()
        nearest_point_line_dict = md.MultiDict()
        for i in range(0, len(source_point_frame)):
            str_source: str = source_point_frame['source'][i].lstrip('(').rstrip(')')
            str_point: str = source_point_frame['point'][i].lstrip('(').rstrip(')')
            source_coord = (float(str_source.split(',')[0]), float(str_source.split(',')[1]))
            point_coord = (float(str_point.split(',')[0]), float(str_point.split(',')[1]))
            source_point_line_dict.put(source_coord, point_coord)
        for i in range(0, len(line_point_frame)):
            line_string_wkt = line_point_frame['nearest_line_wkt'][i]
            str_point: str = line_point_frame['nearest_point'][i].lstrip('(').rstrip(')')
            point_coord = (float(str_point.split(',')[0]), float(str_point.split(',')[1]))
            nearest_point_line_dict.add(line_string_wkt, point_coord)
    else:
        source_point_line_dict = bd.bidict()
        nearest_point_line_dict = md.MultiDict()
        for x in range(0, len(gpd_df_nodes)):
            source_target_nodes_geom = gpd_df_nodes.geometry[x]
            nearest_line = geopandas_min_dist(gpd_df_nodes.geometry[x], gpd_df_network).geometry
            # 用以解决低质数据，下文while处若循环1000次后还连不上，就放弃
            outside_node_grid_geom_to_nearest_p \
                = [q.wkt for q in nearest_points(source_target_nodes_geom, nearest_line)]
            # 从两个最近点中取一个
            nearest_p = (loads(outside_node_grid_geom_to_nearest_p[1]).x, loads(outside_node_grid_geom_to_nearest_p[1]).y)
            source_coord = (source_target_nodes_geom.x, source_target_nodes_geom.y)
            point_src = nearest_line.coords[0]
            point_dst = nearest_line.coords[-1]
            if (nearest_p[0] == point_src[0]) & (nearest_p[1] == point_src[1]):
                point_src_rand = (point_src[0], point_src[1])
                source_point_line_dict.put(source_coord, point_src_rand)
                nearest_point_line_dict.add(nearest_line.wkt, point_src_rand)
            elif (nearest_p[0] == point_dst[0]) & (nearest_p[1] == point_dst[1]):
                point_dst_rand = (point_dst[0], point_dst[1])
                source_point_line_dict.put(source_coord, point_dst_rand)
                nearest_point_line_dict.add(nearest_line.wkt, point_dst_rand)
            else:
                join_line = LineString([source_coord, nearest_p])
                splits_edges_coll: GeometryCollection = geometry.GeometryCollection()
                join_line_length = join_line.length
                if join_line_length == 0:
                    join_line_length = 1e-7
                start_buf_ratio = 1 + (0.001 / join_line_length)
                break_count = 0
                while (len(splits_edges_coll.geoms)) != 2:
                    if join_line.length != 0:
                        join_line_interpolation = get_extrapolated_line(source_coord, nearest_p, start_buf_ratio)
                        splits_edges_coll = split(nearest_line, join_line_interpolation)
                    else:
                        break
                    start_buf_ratio += 0.01 / join_line_length
                    break_count += 1
                    if break_count == 1000:
                        sys.exit('Did over 1000 loops.')
                if join_line.length > 0:
                    nearest_point_exact = splits_edges_coll.geoms[0].coords[-1]
                    nearest_point_exact_rand = (nearest_point_exact[0], nearest_point_exact[1])
                else:
                    nearest_point_exact = nearest_p
                    nearest_point_exact_rand = (nearest_point_exact[0], nearest_point_exact[1])
                source_point_line_dict.put(source_coord, nearest_point_exact_rand)
                nearest_point_line_dict.add(nearest_line.wkt, nearest_point_exact_rand)
        source_point_line_frame = pd.DataFrame({'source': source_point_line_dict.keys(), 'point': source_point_line_dict.values()})
        nearest_point_line_frame = pd.DataFrame({'nearest_line_wkt': nearest_point_line_dict.keys(),
                                                 'nearest_point': nearest_point_line_dict.values()})
        source_point_line_frame.to_csv('source_project_points.csv', index=False)
        nearest_point_line_frame.to_csv('nearest_line_project_points.csv', index=False)
    return source_point_line_dict, nearest_point_line_dict


def build_graph(nodes_df: GeoDataFrame, edges_df: GeoDataFrame):
    """
    参数
    ----------
    nodes_df: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    edges_df: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe

    返回
    ----------
    network_graph: 河网与站点投影点共同构成的图
    """
    if (outdated is False) & (os.path.exists('network_graph.edgelist')):
        network_graph = nx.read_edgelist(
            'network_graph.edgelist', delimiter='|',
            nodetype=lambda node: tuple([float(v) for v in str(node).strip("()").split(",")]),
            create_using=nx.DiGraph)
    else:
        network_graph: DiGraph = nx.DiGraph()
        nearest_point_line_dict = tie_outside_node(nodes_df, edges_df)[1]
        for line in edges_df.geometry:
            line_wkt = line.wkt
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
        nx.write_edgelist(network_graph, 'network_graph.edgelist', delimiter='|')
    return network_graph


def get_upstream_stations(gpd_nodes_df: GeoDataFrame, gpd_network_df: GeoDataFrame,
                          station_index: int, cutoff: int = 2147483647):
    """
    参数
    ----------
    gpd_nodes_df: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    gpd_network_df: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe
    station_index: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以上溯几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量

    返回
    ----------
    upstream_graph: 从当前站点开始，河网与上游所有站点投影点构成的上游子图
    set_up_no_dup: 含有当前站点所有上游站点信息，元素为string形式
    """
    stations_up_list = []
    stations_graph = build_graph(gpd_nodes_df, gpd_network_df)
    source_target_dict = tie_outside_node(gpd_nodes_df, gpd_network_df)[0]
    source_node_coord = (gpd_nodes_df.geometry[station_index].x, gpd_nodes_df.geometry[station_index].y)
    target_node_coord = source_target_dict.get(source_node_coord)
    if os.path.exists("upstream_graph_" + str(station_index) + "_cutoff_" + str(cutoff) + ".edgelist"):
        if (outdated is False) & (os.path.exists('up_down_paths.txt')):
            set_up_no_dup = set()
            with open("up_down_paths.txt", mode='r+') as fp:
                for line_str in fp.readlines():
                    line_str_head = line_str.split(':')[0]
                    line_str_body = line_str.split(':')[1]
                    line_str_stations = line_str_body.split(',')
                    for i in range(0, len(line_str_stations)):
                        line_str_stations[i] = line_str_stations[i].strip('\n').strip(' ').strip("'").strip('"')
                    if (line_str_head == 'upstream') & ('station' + str(station_index) == line_str_stations[-1]):
                        set_up_no_dup.add(str(line_str_stations[-cutoff:]))
            upstream_graph = nx.read_edgelist("upstream_graph_" + str(station_index) + "_cutoff_" + str(cutoff) + ".edgelist",
                                              delimiter='|', nodetype=lambda t: tuple([float(v) for v in str(t).strip("()").split(",")]),
                                              create_using=nx.DiGraph)
        else:
            set_up_no_dup = set()
            upstream_graph = nx.read_edgelist("upstream_graph_" + str(station_index) + "_cutoff_" + str(cutoff) + ".edgelist",
                                              delimiter='|', nodetype=lambda t: tuple([float(v) for v in str(t).strip("()").split(",")]),
                                              create_using=nx.DiGraph)
    else:
        source_target_dict = tie_outside_node(gpd_nodes_df, gpd_network_df)[0]
        source_node_coord = (gpd_nodes_df.geometry[station_index].x, gpd_nodes_df.geometry[station_index].y)
        target_node_coord = source_target_dict.get(source_node_coord)
        upstream_graph = stations_graph.subgraph(nx.ancestors(stations_graph, target_node_coord) | {target_node_coord}).copy()
        if cutoff != 2147483647:
            nx.write_edgelist(upstream_graph,
                              "upstream_graph_" + str(station_index) + "_cutoff_" + str(cutoff) + ".edgelist", delimiter='|')
        set_up_no_dup = set()  # 给找到的路径去重
    for coord in upstream_graph.nodes:
        if upstream_graph.in_degree(coord) == 0 & nx.generic.has_path(upstream_graph, coord, target_node_coord):
            for up_path in nx.all_simple_paths(upstream_graph, coord, target_node_coord):
                for path_coord in up_path:
                    if path_coord in source_target_dict.values():
                        point = Point(source_target_dict.inv[path_coord])
                        for x in range(0, len(gpd_nodes_df)):
                            if gpd_nodes_df.geometry[x] == point:
                                stations_up_list.append('station' + str(x))
                if len(stations_up_list) > 1:
                    set_up_no_dup.add(str(stations_up_list[-cutoff:]))
                stations_up_list.clear()
    return upstream_graph, set_up_no_dup


def get_downstream_stations(gpd_nodes_df: GeoDataFrame, gpd_network_df: GeoDataFrame,
                            station_index: int, cutoff: int = 2147483647):
    """
    参数
    ----------
    gpd_nodes_df: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    gpd_network_df: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe
    station_index: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以往下寻找几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量

    返回
    ----------
    list_stations[:cutoff]: 目标下游站点最近的几个站
    """
    list_stations = []
    if (outdated is False) & os.path.exists('up_down_paths.txt'):
        with open('up_down_paths.txt', mode='r+') as fp:
            for line_str in fp.readlines():
                line_str_head = line_str.split(':')[0]
                line_str_body = line_str.split(':')[1]
                line_str_station = line_str_body.split(',')
                for i in range(0, len(line_str_station)):
                    line_str_station[i] = line_str_station[i].strip('\n').strip(' ').strip("'").strip('"')
                if (line_str_head == 'downstream') & (('station' + str(station_index)) == line_str_station[0]):
                    list_stations = line_str_station
    else:
        stations_graph = build_graph(gpd_nodes_df, gpd_network_df)
        source_target_dict = tie_outside_node(gpd_nodes_df, gpd_network_df)[0]
        source_node_coord = (gpd_nodes_df.geometry[station_index].x, gpd_nodes_df.geometry[station_index].y)
        target_node_coord = source_target_dict.get(source_node_coord)
        dfs_tree_path = nx.dfs_tree(stations_graph, target_node_coord)
        for coord in dfs_tree_path:
            if coord in source_target_dict.values():
                source_coord = source_target_dict.inv[coord]
                point = Point(source_coord)
                for x in range(0, len(gpd_nodes_df)):
                    if gpd_nodes_df.geometry[x] == point:
                        list_stations.append('station' + str(x))
    return list_stations[:cutoff]


def get_upstream_stations_graph(gpd_nodes_df: GeoDataFrame, gpd_network_df: GeoDataFrame, number: int, cutoff: int = 2147483647):
    """
    参数
    ----------
    gpd_nodes_df: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    gpd_network_df: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe
    number: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以上溯几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量

    返回
    ----------
    upstream_graph: 从当前站点开始，河网与上游所有站点投影点构成的上游子图
    new_graph: 当前站点及上游所有站点构成的子图
    origin_graph: 根据河源唯远判断干支流的原始河流抽象图
    """
    upstream_graph = get_upstream_stations(gpd_nodes_df, gpd_network_df, number, cutoff)[0]
    if (outdated is False) & os.path.exists('origin_graph.edgelist'):
        origin_graph = nx.read_edgelist('origin_graph.edgelist', delimiter='|',
                                        nodetype=lambda t: tuple([float(v) for v in str(t).strip("()").split(",")]),
                                        create_using=nx.DiGraph, data=(("weight", float),))
    else:
        origin_graph = nx.DiGraph()
        for line in gpd_network_df.geometry:
            src = line.coords[0]
            dest = line.coords[-1]
            origin_graph.add_edge(src, dest, weight=line.length)
    new_graph = nx.DiGraph()
    source_target_dict = tie_outside_node(gpd_nodes_df, gpd_network_df)[0]
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
    nx.write_edgelist(origin_graph, 'origin_graph.edgelist', delimiter='|', data=["weight"])
    return upstream_graph, new_graph, origin_graph


def write_path_file(gpd_nodes_df: GeoDataFrame, gpd_network_df: GeoDataFrame):
    """
    将上下游信息分行保存在txt文件中以便下次读取，upstream: 代表上游，downstream: 代表下游

    参数
    ----------
    gpd_nodes_df: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    gpd_network_df: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe
    """
    with open("up_down_paths.txt", mode='w+') as fp:
        for i in range(0, len(gpd_nodes_df)):
            print('Writing stream: ' + str(i))
            set_up_no_dup = get_upstream_stations(gpd_nodes_df, gpd_network_df, i)[1]
            list_up_no_dup = list(set_up_no_dup)
            for str_path in list_up_no_dup:
                fp.writelines("upstream: " + str_path.lstrip('[').rstrip(']') + "\n")
            list_down = get_downstream_stations(gpd_nodes_df, gpd_network_df, i)
            fp.writelines("downstream: " + str(list_down).lstrip('[').rstrip(']') + "\n")


def show_upstream_stations_graph(gpd_nodes_df: GeoDataFrame, gpd_network_df: GeoDataFrame, number: int, cutoff: int = 2147483647):
    """
    显示当前站点的上游站点图

    参数
    ----------
    gpd_nodes_df: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    gpd_network_df: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe
    number: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以上溯几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量
    """
    upstream_graph = get_upstream_stations_graph(gpd_nodes_df, gpd_network_df, number, cutoff)[1]
    set_up_no_dup = get_upstream_stations(gpd_nodes_df, gpd_network_df, number, cutoff)[1]
    for list_str in set_up_no_dup:
        print(list_str.lstrip('[').rstrip(']'))  # 输出的都是字串
    nx.draw(upstream_graph, node_size=10)
    plt.show()


def show_downstream_stations(gpd_nodes_df: GeoDataFrame, gpd_network_df: GeoDataFrame, number: int, cutoff: int = 2147483647):
    """
    显示当前站点的下游站点列表

    参数
    ----------
    gpd_nodes_df: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    gpd_network_df: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe
    number: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以往下寻找几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量
    """
    list_stations = get_downstream_stations(gpd_nodes_df, gpd_network_df, number, cutoff)
    print(list_stations)


def upstream_node_on_mainstream(gpd_nodes_df, gpd_network_df, number_src, number_target):
    """
    以number_src为当前站点，判断number_target所代表的站点是否存在于当前站点上游流域的干支流中，存在三种结果：
    在支流中（In tributary），在干流中（In Mainstream），不在上游流域中

    参数
    ----------
    gpd_nodes_df: 存储站点数据的GeoDataFrame，参考本文件开头gpd_nodes_dataframe
    gpd_network_df: 存储河网（所有LineString）数据的GeoDataFrame，参考本文件开头gpd_network_dataframe
    number_src: 用以确定上游流域的站点号，可以从图层属性表中辅助判断
    number_target: 待判断的站点号
    """
    # number_src是要生成子图的原点号，number_target是要判断干支流的点号
    source_point = (gpd_nodes_df.geometry[number_src].x, gpd_nodes_df.geometry[number_src].y)
    target_point = (gpd_nodes_df.geometry[number_target].x, gpd_nodes_df.geometry[number_target].y)
    if (outdated is False) & (os.path.exists('nearest_line_project_points.csv') & os.path.exists('source_project_points.csv')):
        nearest_line_project_df = pd.read_csv('nearest_line_project_points.csv')
        source_project_point_dict = tie_outside_node(gpd_nodes_df, gpd_network_df)[0]
        source_nearest_point = source_project_point_dict.get(source_point)
        target_nearest_point = source_project_point_dict.get(target_point)
        nearest_source_line_wkt = nearest_line_project_df['nearest_line_wkt'][nearest_line_project_df['nearest_point'] ==
                                                                              str(source_nearest_point)].values[0]
        nearest_target_line_wkt = nearest_line_project_df['nearest_line_wkt'][nearest_line_project_df['nearest_point'] ==
                                                                              str(target_nearest_point)].values[0]
        nearest_source_line: LineString = loads(nearest_source_line_wkt)
        nearest_target_line: LineString = loads(nearest_target_line_wkt)
    else:
        nearest_target_line: LineString = geopandas_min_dist(Point(target_point), gpd_network_df).geometry
        nearest_source_line: LineString = geopandas_min_dist(Point(source_point), gpd_network_df).geometry
    origin_graph = get_upstream_stations_graph(gpd_nodes_df, gpd_network_df, number_src)[2]
    origin_target_point = nearest_target_line.coords[0]
    origin_source_point = nearest_source_line.coords[-1]
    origin_graph_src_sub: DiGraph = nx.subgraph(origin_graph, nx.ancestors(origin_graph, origin_source_point) | {
        origin_source_point}).copy()
    set_up_no_dup = get_upstream_stations(gpd_nodes_df, gpd_network_df, number_src, cutoff=10000)[1]
    in_basin = False
    for upstream_str in set_up_no_dup:
        in_basin = in_basin | ('station' + str(number_target) in upstream_str)
        if in_basin is False:
            continue
        else:
            # 寻找所有到nearest_line.coords[-1]的边，然后判断其是否小于DAG最长长度
            line_stations = upstream_str.strip('[]').replace("'", "").replace(' ', '').split(',')
            for i in range(0, len(line_stations)):
                line_stations[i] = line_stations[i].replace('"', '')
            if (nx.shortest_path_length(origin_graph_src_sub, origin_target_point, origin_source_point, weight='weight') <
                nx.dag_longest_path_length(origin_graph_src_sub, weight='weight')) & (origin_target_point not in nx.dag_longest_path(
                origin_graph_src_sub, weight='weight')):
                print('In tributary: ' + str(line_stations))
                break
            else:
                cutoff = line_stations.index('station' + str(number_target))
                print('In Mainstream: ' + str(line_stations[cutoff:]))
                break
    if in_basin is False:
        print(str(number_target) + ' is not in upstream basin of ' + str(number_src))


if __name__ == '__main__':
    index = 4
    START_TIME = time.time()
    show_upstream_stations_graph(gpd_nodes_dataframe, gpd_network_dataframe, index, 5)
    print('____________________________________________________________________________')
    print(upstream_node_on_mainstream(gpd_nodes_dataframe, gpd_network_dataframe, index, 28))
    print('____________________________________________________________________________')
    show_downstream_stations(gpd_nodes_dataframe, gpd_network_dataframe, index)
    if outdated is True:
        write_path_file(gpd_nodes_dataframe, gpd_network_dataframe)
    STOP_TIME = time.time()
    TOTAL_TIME = STOP_TIME - START_TIME
    print('done:', TOTAL_TIME, 'seconds')
    print("Finished")
