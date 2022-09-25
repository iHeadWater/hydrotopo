import os
import sys

import bidict as bd
import click
import matplotlib.pyplot as plt
import multidict as md
import networkx as nx
import pandas as pd
from networkx import DiGraph
from shapely import geometry
from shapely.geometry import Point, LineString, GeometryCollection
from shapely.ops import nearest_points, split
from shapely.wkt import loads
from shapefile import Reader, Shape

"""
本文件线元测试数据来自https://www.hydrosheds.org/products/hydrorivers，坐标为EPSG:4326，单位全部为度，1°≈110km
河网矢量图层应在使用前在GIS软件中转换为单部件，否则可能出现问题
"""


def shape_type_to_wkt(shape: Shape):
    geo_shpfile = shape.__geo_interface__
    geo_class = geo_shpfile['type'].upper()
    wkt: str = geo_class + ' ('
    geo_coords = []
    if geo_class == 'POINT':
        geo_coords.append(geo_shpfile['coordinates'])
    else:
        geo_coords = list(geo_shpfile['coordinates'])
    for coord_tuple in geo_coords:
        wkt += str(coord_tuple[0])
        wkt += ' '
        wkt += str(coord_tuple[1])
        wkt += ', '
    wkt = wkt.rstrip(', ')
    wkt += ')'
    return wkt


def line_min_dist(pyshp_coord, net_reader: Reader, initial_buffer=0.005):
    """
    https://gis.stackexchange.com/questions/266730/filter-by-bounding-box-in-geopandas/266833
    从限定范围内获取图元，本例中为河网线

    参数
    ----------
    point : tuple
        待寻找最近线段的点
    net_reader : shapefile.Reader
        存储河网线段的Reader数据结构
    initial_buffer : float
        给点加的buffer（搜寻半径）大小，初始为0.005度，约为550米

    返回
    ----------
    离point最近的一条Shapely LineString
    """
    shapely_point = Point(pyshp_coord)
    buffer_steps = 0.002
    # 给点或者其他geometry加buffer
    bbox = shapely_point.buffer(initial_buffer).bounds
    shapes = net_reader.shapes(bbox)
    while len(shapes) == 0:
        initial_buffer = initial_buffer + buffer_steps
        bbox = (shapely_point.buffer(initial_buffer)).bounds
        shapes = net_reader.shapes(bbox)
    min_dist = 360
    min_dis_line = shapes[0]
    for shape in shapes:
        shapely_shape = loads(shape_type_to_wkt(shape))
        if shapely_point.distance(shapely_shape) < min_dist:
            min_dist = shapely_point.distance(shapely_shape)
            min_dis_line = shape
    return loads(shape_type_to_wkt(min_dis_line))


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


def tie_outside_node(node_reader: Reader, network_reader: Reader, outdated: bool, cache_dir=os.curdir, output_path=os.curdir):
    """
    值得注意，由于tie_outside_node操作很费时，故工程中如有source_project_points.csv文件和nearest_line_project_points.csv文件，
    程序会优先从这两个文件中读取数据，所以如要从头生成，必须删除这两个文件

    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader

    返回
    ----------
    nearest_point_line_dict: 线上投影点和最近河流线的对应列表
    source_point_line_dict: 原站点和站点在线上投影点的对应列表
    """
    if (outdated is False) | ((outdated is True) & os.path.exists(os.path.join(cache_dir, 'source_project_points.csv')) &
                              os.path.exists(os.path.join(cache_dir, 'nearest_line_project_points.csv'))):
        source_point_frame = pd.read_csv(os.path.join(cache_dir, 'source_project_points.csv'))
        line_point_frame = pd.read_csv(os.path.join(cache_dir, 'nearest_line_project_points.csv'))
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
        for x in range(0, len(node_reader)):
            source_coord = tuple(node_reader.shape(x).points[0])
            nearest_line: LineString = line_min_dist(source_coord, network_reader)
            # 用以解决低质数据，下文while处若循环1000次后还连不上，就放弃
            outside_node_grid_geom_to_nearest_p \
                = [q.wkt for q in nearest_points(Point(source_coord), nearest_line)]
            # 从两个最近点中取一个
            nearest_p = (loads(outside_node_grid_geom_to_nearest_p[1]).x, loads(outside_node_grid_geom_to_nearest_p[1]).y)
            point_src = nearest_line.coords[0]
            point_dst = nearest_line.coords[-1]
            if (nearest_p[0] == point_src[0]) & (nearest_p[1] == point_src[1]):
                gradient = (point_dst[1] - point_src[1]) / (point_dst[0] - point_src[0])
                point_src_rand = (point_src[0] + 1e-6, point_src[1] + gradient * 1e-6)
                source_point_line_dict.put(source_coord, point_src_rand)
                nearest_point_line_dict.add(nearest_line.wkt, point_src_rand)
            elif (nearest_p[0] == point_dst[0]) & (nearest_p[1] == point_dst[1]):
                gradient = (point_dst[1] - point_src[1]) / (point_dst[0] - point_src[0])
                point_dst_rand = (point_dst[0] + 1e-6, point_dst[1] + gradient * 1e-6)
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
                    nearest_point_exact = (nearest_point_exact[0], nearest_point_exact[1])
                else:
                    nearest_point_exact = nearest_p
                    nearest_point_exact = (nearest_point_exact[0], nearest_point_exact[1])
                source_point_line_dict.put(source_coord, nearest_point_exact)
                nearest_point_line_dict.add(nearest_line.wkt, nearest_point_exact)
        source_point_line_frame = pd.DataFrame({'source': source_point_line_dict.keys(), 'point': source_point_line_dict.values()})
        nearest_point_line_frame = pd.DataFrame({'nearest_line_wkt': nearest_point_line_dict.keys(),
                                                 'nearest_point': nearest_point_line_dict.values()})
        source_point_line_frame.to_csv(os.path.join(output_path, 'source_project_points.csv'), index=False)
        nearest_point_line_frame.to_csv(os.path.join(output_path, 'nearest_line_project_points.csv'), index=False)
    return source_point_line_dict, nearest_point_line_dict


def build_graph(node_reader: Reader, network_reader: Reader, outdated: bool, cache_dir=os.curdir, output_path=os.curdir):
    """
    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader

    返回
    ----------
    network_graph: 河网与站点投影点共同构成的图
    """
    if (outdated is False) & (os.path.exists(os.path.join(cache_dir, 'network_graph.edgelist'))):
        network_graph = nx.read_edgelist(
            os.path.join(cache_dir, 'network_graph.edgelist'), delimiter='|',
            nodetype=lambda node: tuple([float(v) for v in str(node).strip("()").split(",")]),
            create_using=nx.DiGraph)
    else:
        network_graph: DiGraph = nx.DiGraph()
        nearest_point_line_dict = tie_outside_node(node_reader, network_reader, outdated, cache_dir, output_path)[1]
        for line in network_reader.shapes():
            line_wkt: str = shape_type_to_wkt(line)
            list_point_and_dis = []
            list_coords = []
            if line_wkt in nearest_point_line_dict.keys():
                for coord in nearest_point_line_dict.getall(line_wkt):
                    list_point_and_dis.append((coord, Point(coord).distance(Point(line.points[0]))))
                list_point_and_dis.append((line.points[0], 0))
                list_point_and_dis.append((line.points[-1], 360))
                list_point_and_dis = sorted(list_point_and_dis, key=lambda point_dis: point_dis[1])
                for i in range(0, len(list_point_and_dis)):
                    list_coords.append(list_point_and_dis[i][0])
                nx.add_path(network_graph, list_coords)
            else:
                src = line.points[0]
                dest = line.points[-1]
                network_graph.add_edge(src, dest)
        nx.write_edgelist(network_graph, os.path.join(output_path, 'network_graph.edgelist'), delimiter='|')
    return network_graph


def get_upstream_stations(node_reader: Reader, network_reader: Reader,
                          station_index: int, outdated: bool, cutoff: int = 2147483647, cache_dir=os.curdir, output_path=os.curdir):
    """
    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader
    station_index: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以上溯几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量

    返回
    ----------
    upstream_graph: 从当前站点开始，河网与上游所有站点投影点构成的上游子图
    set_up_no_dup: 含有当前站点所有上游站点信息，元素为string形式
    """
    stations_up_list = []
    stations_graph = build_graph(node_reader, network_reader, outdated, cache_dir, output_path)
    source_target_dict = tie_outside_node(node_reader, network_reader, outdated, cache_dir, output_path)[0]
    source_node_coord = tuple(node_reader.shape(station_index).points[0])
    target_node_coord = source_target_dict.get(source_node_coord)
    upstream_graph_path = "upstream_graph_" + str(station_index) + "_cutoff_" + str(cutoff) + ".edgelist"
    if os.path.exists(os.path.join(cache_dir, upstream_graph_path)):
        if (outdated is False) & (os.path.exists(os.path.join(cache_dir, 'up_down_paths.txt'))):
            set_up_no_dup = set()
            with open(os.path.join(cache_dir, 'up_down_paths.txt'), mode='r+') as fp:
                for line_str in fp.readlines():
                    line_str_head = line_str.split(':')[0]
                    line_str_body = line_str.split(':')[1]
                    line_str_stations = line_str_body.split(',')
                    for i in range(0, len(line_str_stations)):
                        line_str_stations[i] = line_str_stations[i].strip('\n').strip(' ').strip("'").strip('"')
                    if (line_str_head == 'upstream') & ('sta' + str(station_index) == line_str_stations[-1]):
                        set_up_no_dup.add(str(line_str_stations[-cutoff:]))
            upstream_graph = nx.read_edgelist(os.path.join(cache_dir, upstream_graph_path), delimiter='|',
                                              nodetype=lambda t: tuple([float(v) for v in str(t).strip("()").split(",")]),
                                              create_using=nx.DiGraph)
        else:
            set_up_no_dup = set()
            upstream_graph = nx.read_edgelist(os.path.join(cache_dir, upstream_graph_path), delimiter='|',
                                              nodetype=lambda t: tuple([float(v) for v in str(t).strip("()").split(",")]),
                                              create_using=nx.DiGraph)
    else:
        source_target_dict = tie_outside_node(node_reader, network_reader, outdated, cache_dir, output_path)[0]
        source_node_coord = tuple(node_reader.shape(station_index).points[0])
        target_node_coord = source_target_dict.get(source_node_coord)
        upstream_graph = stations_graph.subgraph(nx.ancestors(stations_graph, target_node_coord) | {target_node_coord}).copy()
        if cutoff != 2147483647:
            nx.write_edgelist(upstream_graph,
                              os.path.join(output_path, upstream_graph_path), delimiter='|')
        set_up_no_dup = set()  # 给找到的路径去重
    for coord in upstream_graph.nodes:
        if upstream_graph.in_degree(coord) == 0 & nx.generic.has_path(upstream_graph, coord, target_node_coord):
            for up_path in nx.all_simple_paths(upstream_graph, coord, target_node_coord):
                for path_coord in up_path:
                    if path_coord in source_target_dict.values():
                        coord = source_target_dict.inv[path_coord]
                        for x in range(0, len(node_reader)):
                            if tuple(node_reader.shape(x).points[0]) == coord:
                                stations_up_list.append('sta' + str(x))
                if len(stations_up_list) > 1:
                    set_up_no_dup.add(str(stations_up_list[-cutoff:]))
                stations_up_list.clear()
    return upstream_graph, set_up_no_dup


def get_downstream_stations(node_reader: Reader, network_reader: Reader,
                            station_index: int, outdated: bool, cutoff: int = 2147483647, cache_dir=os.curdir, output_path=os.curdir):
    """
    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader
    station_index: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以往下寻找几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量

    返回
    ----------
    list_stations[:cutoff]: 目标下游站点最近的几个站
    """
    list_stations = []
    if (outdated is False) & os.path.exists(os.path.exists(os.path.join(cache_dir, 'up_down_paths.txt'))):
        with open(os.path.join(cache_dir, 'up_down_paths.txt'), mode='r+') as fp:
            for line_str in fp.readlines():
                line_str_head = line_str.split(':')[0]
                line_str_body = line_str.split(':')[1]
                line_str_station = line_str_body.split(',')
                for i in range(0, len(line_str_station)):
                    line_str_station[i] = line_str_station[i].strip('\n').strip(' ').strip("'").strip('"')
                if (line_str_head == 'downstream') & (('sta' + str(station_index)) == line_str_station[0]):
                    list_stations = line_str_station
    else:
        stations_graph = build_graph(node_reader, network_reader, outdated, cache_dir, output_path)
        source_target_dict = tie_outside_node(node_reader, network_reader, outdated, cache_dir, output_path)[0]
        source_node_coord = tuple(node_reader.shape(station_index).points[0])
        target_node_coord = source_target_dict.get(source_node_coord)
        dfs_tree_path = nx.dfs_tree(stations_graph, target_node_coord)
        for coord in dfs_tree_path:
            if coord in source_target_dict.values():
                source_coord = source_target_dict.inv[coord]
                for x in range(0, len(node_reader)):
                    if tuple(node_reader.shape(x).points[0]) == source_coord:
                        list_stations.append('sta' + str(x))
    return list_stations[:cutoff]


def get_upstream_stations_graph(node_reader: Reader, network_reader: Reader, number: int, outdated: bool,
                                cutoff: int = 2147483647, cache_dir=os.curdir, output_path=os.curdir):
    """
    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader
    number: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以上溯几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量

    返回
    ----------
    upstream_graph: 从当前站点开始，河网与上游所有站点投影点构成的上游子图
    new_graph: 当前站点及上游所有站点构成的子图
    origin_graph: 根据河源唯远判断干支流的原始河流抽象图
    """
    upstream_graph = get_upstream_stations(node_reader, network_reader, number, outdated, cutoff, cache_dir, output_path)[0]
    if (outdated is False) & os.path.exists(os.path.join(cache_dir, 'origin_graph.edgelist')):
        origin_graph = nx.read_edgelist(os.path.join(cache_dir, 'origin_graph.edgelist'), delimiter='|',
                                        nodetype=lambda t: tuple([float(v) for v in str(t).strip("()").split(",")]),
                                        create_using=nx.DiGraph, data=(("weight", float),))
    else:
        origin_graph = nx.DiGraph()
        for line in network_reader.shapes():
            src = tuple(line.points[0])
            dest = tuple(line.points[-1])
            shapely_line: LineString = loads(shape_type_to_wkt(line))
            origin_graph.add_edge(src, dest, weight=shapely_line.length)
    new_graph = nx.DiGraph()
    source_target_dict = tie_outside_node(node_reader, network_reader, outdated, cache_dir, output_path)[0]
    source_node_coord = tuple(node_reader.shape(number).points[0])
    target_node_coord = source_target_dict.get(source_node_coord)
    for coord in upstream_graph.nodes:
        if upstream_graph.in_degree(coord) == 0:
            for up_path in nx.all_simple_paths(upstream_graph, coord, target_node_coord):
                list_up_path = list(up_path)
                for path_coord in up_path:
                    if path_coord in origin_graph.nodes:
                        list_up_path.remove(path_coord)
                nx.add_path(new_graph, list_up_path[-cutoff:])
    nx.write_edgelist(origin_graph, os.path.join(output_path, 'origin_graph.edgelist'), delimiter='|', data=["weight"])
    return upstream_graph, new_graph, origin_graph


def write_path_file(node_reader: Reader, network_reader: Reader, output_path=os.curdir):
    """
    将上下游信息分行保存在txt文件中以便下次读取，upstream: 代表上游，downstream: 代表下游

    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader
    """
    with open(os.path.join(output_path, "up_down_paths.txt"), mode='w+') as fp:
        for i in range(0, len(node_reader)):
            click.echo('Writing stream: ' + str(i))
            set_up_no_dup = get_upstream_stations(node_reader, network_reader, i, True, 2147483647, output_path)[1]
            list_up_no_dup = list(set_up_no_dup)
            for str_path in list_up_no_dup:
                fp.writelines("upstream: " + str_path.lstrip('[').rstrip(']') + "\n")
            list_down = get_downstream_stations(node_reader, network_reader, i, True)
            fp.writelines("downstream: " + str(list_down).lstrip('[').rstrip(']') + "\n")


def show_upstream_stations_graph(node_reader: Reader, network_reader: Reader, number: int, outdated: bool,
                                 cache_dir=os.curdir, output_path=os.curdir, cutoff: int = 2147483647):
    """
    显示当前站点的上游站点图

    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader
    number: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以上溯几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量
    """
    upstream_graph = get_upstream_stations_graph(node_reader, network_reader, number, outdated, cutoff, cache_dir, output_path)[1]
    set_up_no_dup = get_upstream_stations(node_reader, network_reader, number, outdated, cutoff, cache_dir, output_path)[1]
    for list_str in set_up_no_dup:
        click.echo(list_str.lstrip('[').rstrip(']'))  # 输出的都是字串
    nx.draw_networkx(upstream_graph, node_size=10)
    plt.savefig(os.path.join(output_path, 'upstream_graph_' + str(number) + '_cutoff_' + str(cutoff) + '.png'))
    plt.show()


def show_downstream_stations(node_reader: Reader, network_reader: Reader, number: int, outdated: bool,
                             cache_dir=os.curdir, output_path=os.curdir, cutoff: int = 2147483647):
    """
    显示当前站点的下游站点列表

    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader
    number: 站点图层中的标号，本例中从0开始，针对不同数据，可以从图层属性表中行号辅助判断
    cutoff: 同一条河上最多可以往下寻找几个水文站，默认为int.max（当然也可以指定一个很大的数），即不限数量
    """
    list_stations = get_downstream_stations(node_reader, network_reader, number, outdated, cutoff, cache_dir, output_path)
    click.echo(list_stations)


def upstream_node_on_mainstream(node_reader: Reader, network_reader: Reader, number_src, number_target, outdated: bool,
                                cache_dir=os.curdir):
    """
    以number_src为当前站点，判断number_target所代表的站点是否存在于当前站点上游流域的干支流中，存在三种结果：
    在支流中（In tributary），在干流中（In Mainstream），不在上游流域中

    参数
    ----------
    node_reader: 存储站点数据的Reader
    network_reader: 存储河网（所有LineString）数据的Reader
    number_src: 用以确定上游流域的站点号，可以从图层属性表中辅助判断
    number_target: 待判断的站点号
    """
    # number_src是要生成子图的原点号，number_target是要判断干支流的点号
    source_point = tuple(node_reader.shape(int(number_src)).points[0])
    target_point = tuple(node_reader.shape(int(number_target)).points[0])
    if (outdated is False) & (os.path.exists(os.path.join(cache_dir, 'nearest_line_project_points.csv')) & os.path.exists(
            os.path.join(cache_dir, 'source_project_points.csv'))):
        nearest_line_project_df = pd.read_csv(os.path.join(cache_dir, 'nearest_line_project_points.csv'))
        source_project_point_dict = tie_outside_node(node_reader, network_reader, outdated, cache_dir)[0]
        source_nearest_point = source_project_point_dict.get(source_point)
        target_nearest_point = source_project_point_dict.get(target_point)
        nearest_source_line_wkt = nearest_line_project_df['nearest_line_wkt'][nearest_line_project_df['nearest_point'] ==
                                                                              str(source_nearest_point)].values[0]
        nearest_target_line_wkt = nearest_line_project_df['nearest_line_wkt'][nearest_line_project_df['nearest_point'] ==
                                                                              str(target_nearest_point)].values[0]
        nearest_source_line: LineString = loads(nearest_source_line_wkt)
        nearest_target_line: LineString = loads(nearest_target_line_wkt)
    else:
        nearest_target_line: LineString = line_min_dist(target_point, network_reader)
        nearest_source_line: LineString = line_min_dist(source_point, network_reader)
    origin_graph = get_upstream_stations_graph(node_reader, network_reader, number_src, outdated)[2]
    origin_target_point = nearest_target_line.coords[0]
    origin_source_point = nearest_source_line.coords[-1]
    origin_graph_src_sub: DiGraph = nx.subgraph(origin_graph, nx.ancestors(origin_graph, origin_source_point) | {
        origin_source_point}).copy()
    set_up_no_dup = get_upstream_stations(node_reader, network_reader, number_src, outdated, cutoff=10000)[1]
    in_basin = False
    for upstream_str in set_up_no_dup:
        in_basin = in_basin | ('sta' + str(number_target) in upstream_str)
        if in_basin is False:
            continue
        else:
            # 寻找所有到nearest_line.coords[-1]的边，然后判断其是否小于DAG最长长度
            line_stations = upstream_str.strip('[]').replace("'", "").replace(' ', '').split(',')
            for i in range(0, len(line_stations)):
                line_stations[i] = line_stations[i].replace('"', '')
            if (nx.shortest_path_length(origin_graph_src_sub, origin_target_point, origin_source_point, weight='weight') <
                nx.dag_longest_path_length(origin_graph_src_sub, weight='weight')) & (
                    origin_target_point not in nx.dag_longest_path(origin_graph_src_sub, weight='weight')):
                click.echo('In tributary: ' + str(line_stations))
                break
            else:
                cutoff = line_stations.index('sta' + str(number_target))
                click.echo('In Mainstream: ' + str(line_stations[cutoff:]))
                break
    if in_basin is False:
        click.echo(str(number_target) + ' is not in upstream basin of ' + str(number_src))
