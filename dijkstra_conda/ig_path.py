from os.path import relpath

import geopandas as gpd
import igraph as ig
import numpy as np
import shapely
from geopandas import GeoDataFrame
from numpy import ndarray
from shapely import STRtree, MultiLineString, intersection
from shapely.geometry import Point, LineString
from shapely.ops import split, nearest_points

gpd_df_node = gpd.read_file(relpath('test_data/near_changchun_dots.shp'), engine='pyogrio')
gpd_df_network = gpd.read_file(relpath('test_data/near_changchun_cut.shp'), engine='pyogrio')


def get_extrapolated_line(source_coord, coord, extrapolate_ratio):
    source = Point(source_coord)
    p2 = Point(coord)
    c = (source.x + extrapolate_ratio * (p2.x - source.x), source.y + extrapolate_ratio * (p2.y - source.y))
    return LineString([source, c])


def line_min_dist(gpd_node_df: GeoDataFrame, gpd_line_df: GeoDataFrame):
    gpd_min_line_df = gpd.sjoin_nearest(gpd_node_df, gpd_line_df, 'left')
    gpd_min_line_df = gpd_min_line_df.loc[~gpd_min_line_df.index.duplicated()]
    index_right = gpd_min_line_df['index_right'].to_numpy()
    origin_geom_array = gpd_line_df.geometry.to_numpy()
    geom_array = np.delete(origin_geom_array, index_right)
    index_geom_array = np.empty(len(gpd_node_df), dtype=LineString)
    new_geom_array = np.array([])
    past_indexes = np.zeros(shape=len(gpd_min_line_df), dtype=int)
    past_indexes[past_indexes == 0] = -1
    node_indexes = np.arange(0, len(gpd_min_line_df), 1)
    for index in node_indexes:
        if index in past_indexes:
            continue
        else:
            right = index_right[index]
            point_ids = list(gpd_min_line_df.index[gpd_min_line_df['index_right'] == right])
            origin_line = origin_geom_array[right]
            ex_list = []
            inter_dict = {}
            for point_id in point_ids:
                point = gpd_min_line_df.geometry[point_id]
                pre_near_point = sure_nearest_point(point, origin_line)
                line = LineString([point, pre_near_point])
                ex_line = get_extrapolated_line(point, pre_near_point, 10 / line.length)
                inter_point = intersection(ex_line, origin_line)
                inter_dict[point_id] = inter_point
                ex_list.append(ex_line)
                past_indexes[index] = point_id
            multiline = MultiLineString(ex_list)
            min_lines = split(origin_line, multiline)
            geom_list = list(min_lines.geoms)
            new_geom_array = np.append(new_geom_array, geom_list)
            for key in inter_dict.keys():
                for cursor in range(1, min(len(geom_list), len(point_ids) + 1)):
                    if inter_dict[key].geom_type == 'Point':
                        checkpoint = inter_dict[key].coords[0]
                    else:
                        coord_list = list(inter_dict[key].geoms)
                        checkpoint = coord_list[0].coords[0]
                    if geom_list[cursor].coords[0] == checkpoint:
                        index_geom_array[key] = geom_list[cursor]
    # new_geom_array和index_geom_array不一样长
    geom_array = np.append(geom_array, new_geom_array)
    return geom_array, new_geom_array, index_geom_array


def sure_nearest_point(source_coord: Point, origin_line: LineString):
    nearest_p = nearest_points(source_coord, origin_line)[1]
    point_src = origin_line.coords[0]
    point_dst = origin_line.coords[-1]
    point_neigh = origin_line.coords[1]
    point_last = origin_line.coords[-2]
    if (nearest_p.x == point_src[0]) & (nearest_p.y == point_src[1]):
        nearest_point = Point((point_src[0] + point_neigh[0]) / 2, (point_src[1] + point_neigh[1]) / 2)
    elif (nearest_p.x == point_dst[0]) & (nearest_p.y == point_dst[1]):
        nearest_point = Point((point_dst[0] + point_last[0]) / 2, (point_src[1] + point_last[1]) / 2)
    elif (source_coord.x == nearest_p.x) & (source_coord.y == nearest_p.y):
        # 如果源点正好在线上，那么投影点就可以是它自身，这时最近点选任何一个点外线都可以，反正交点还是源点
        nearest_point = Point((point_src[0] + point_dst[0]) / 2 + 0.001, (point_src[1] + point_dst[1]) / 2 + 0.001)
    else:
        nearest_point = shapely.line_interpolate_point(origin_line,
                                                       shapely.line_locate_point(origin_line, source_coord))
    return nearest_point


def build_graph(geom_array: ndarray):
    geom_tree = STRtree(geom_array)
    left, right = geom_tree.query(geom_array, predicate='intersects')
    pairs = np.array([left[left != right], right[left != right]]).T
    directed_pairs = np.zeros(shape=2 * len(pairs), dtype=int)
    directed_pairs[directed_pairs == 0] = -1
    # 无向图要修正成有向图
    for i in np.arange(0, len(pairs), 1):
        end0 = geom_array[pairs[i][0]].coords[-1]
        start1 = geom_array[pairs[i][1]].coords[0]
        if end0 == start1:
            directed_pairs[2 * i] = pairs[i][0]
            directed_pairs[2 * i + 1] = pairs[i][1]
    directed_pairs = directed_pairs[directed_pairs != -1]
    directed_pairs = directed_pairs.reshape(-1, 2)
    total_graph = ig.Graph(directed_pairs, directed=True)
    return total_graph


def find_edge_nodes(gpd_nodes_df, gpd_network_df, station_index: int, switch='up', cutoff: int = 2147483647):
    geom_array, new_geom_array, index_geom_array = line_min_dist(gpd_nodes_df, gpd_network_df)
    graph = build_graph(geom_array)
    # 当前站点所对应线索引
    cur_index = np.argwhere(new_geom_array == index_geom_array[station_index])[0][0]
    true_index = len(geom_array) - len(new_geom_array) + cur_index
    # paths里面的是元素所在的（顶点）标号，而非元素本身
    if switch == 'up':
        paths = graph.get_all_shortest_paths(v=true_index, mode='in')
        sta_lists = []
        for path in paths:
            sta_list = []
            for line in path:
                sta_index = np.argwhere(index_geom_array == geom_array[line])
                if len(sta_index) > 0:
                    sta_list.append(sta_index[0][0])
            sta_list.reverse()
            sta_lists.append(sta_list[-cutoff:])
    elif switch == 'down':
        paths = graph.get_all_shortest_paths(v=true_index, mode='out')
        sta_lists = []
        for path in paths:
            sta_list = []
            for line in path:
                sta_index = np.argwhere(index_geom_array == geom_array[line])
                if len(sta_index) > 0:
                    sta_list.append(sta_index[0][0])
            sta_lists.append(sta_list[:cutoff])
    else:
        sta_lists = []
    return np.unique(np.array(sta_lists))


def calc_distance(gpd_nodes_df, gpd_network_df, start: int, end: int):
    geom_array, new_geom_array, index_geom_array = line_min_dist(gpd_nodes_df, gpd_network_df)
    graph = build_graph(geom_array)
    # 当前站点所对应线索引
    cur_start_index = np.argwhere(new_geom_array == index_geom_array[start])[0][0]
    cur_end_index = np.argwhere(new_geom_array == index_geom_array[end])[0][0]
    true_start_index = len(geom_array) - len(new_geom_array) + cur_start_index
    true_end_index = len(geom_array) - len(new_geom_array) + cur_end_index
    paths = graph.get_all_shortest_paths(v=true_start_index, to=true_end_index, mode='all')
    len_list = []
    for path in paths:
        total_len = shapely.distance(gpd_nodes_df.geometry[start], geom_array[path[0]])
        for number in path:
            total_len += (geom_array[number]).length
        total_len += shapely.distance(gpd_nodes_df.geometry[end], geom_array[path[-1]])
        len_list.append(total_len)
    return len_list


def test_upper_nodes():
    print(find_edge_nodes(gpd_df_node, gpd_df_network, 22, 'up', 6))
    print('__________________________________________________')
    print(find_edge_nodes(gpd_df_node, gpd_df_network, 10, 'down', 6))


def test_calc_distance():
    # 倒过来也可以成立
    print(calc_distance(gpd_df_node, gpd_df_network, 10, 7))


def test_down_array():
    print(line_min_dist(gpd_df_node, gpd_df_network)[1])


def test_pairs():
    print(build_graph(line_min_dist(gpd_df_node, gpd_df_network)[0]))
