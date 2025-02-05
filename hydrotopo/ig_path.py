import geopandas as gpd
import igraph as ig
import numpy as np
import shapely
from geopandas import GeoDataFrame
from shapely import STRtree, MultiLineString, intersection
from shapely.geometry import Point, LineString
from shapely.ops import split, nearest_points


def get_extrapolated_line(source_coord, coord, extrapolate_ratio):
    source = Point(source_coord)
    p2 = Point(coord)
    c = (source.x + extrapolate_ratio * (p2.x - source.x), source.y + extrapolate_ratio * (p2.y - source.y))
    return LineString([source, c])


def line_min_dist(gpd_node_df: GeoDataFrame, gpd_line_df: GeoDataFrame):
    # 根据gpd_node_df，确定gpd_line_df中离站点最近的河网线
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
        if index not in past_indexes:
            right = index_right[index]
            # 找出同一条河网线附近一条上下游上的所有站点
            point_ids = list(gpd_min_line_df.index[gpd_min_line_df['index_right'] == right])
            origin_line = origin_geom_array[right]
            ex_list, inter_dict = gen_split_lines(gpd_min_line_df, point_ids, origin_line, past_indexes, index)
            # 将多条插分线组成一个MultiLineString，一次将河网线分开
            multiline = MultiLineString(ex_list)
            min_lines = split(origin_line, multiline)
            geom_list = list(min_lines.geoms)
            new_geom_array = np.append(new_geom_array, geom_list)
            # 由于new_geom_array和index_geom_array不一样长，所以要使用一个index_geom_array记录点索引对应的线索引
            sure_index_geom_array(index_geom_array, inter_dict, geom_list)
    geom_array = np.append(geom_array, new_geom_array)
    return geom_array, new_geom_array, index_geom_array

def sure_index_geom_array(index_geom_array, inter_dict, geom_list):
    for key in inter_dict.keys():
        for cursor in range(1, len(geom_list)):
            if (inter_dict[key].geom_type == 'Point') | (inter_dict[key].geom_type == 'LineString'):
                checkpoint = inter_dict[key].coords[0]
            else:
                coord_list = list(inter_dict[key].geoms)
                checkpoint = coord_list[0].coords[0]
            if geom_list[cursor].coords[0] == checkpoint:
                index_geom_array[key] = geom_list[cursor]
                break

def gen_split_lines(gpd_min_line_df, point_ids, origin_line, past_indexes, index):
    ex_list = []
    inter_dict = {}
    # 针对多个点生成多条插分线
    for point_id in point_ids:
        point = gpd_min_line_df.geometry[point_id]
        pre_near_point = sure_nearest_point(point, origin_line)
        line = LineString([point, pre_near_point])
        ex_line = get_extrapolated_line(point, pre_near_point, 10 / line.length)
        inter_point = intersection(ex_line, origin_line)
        inter_dict[point_id] = inter_point
        ex_list.append(ex_line)
        past_indexes[index] = point_id
    return ex_list, inter_dict

def sure_nearest_point(source_coord: Point, origin_line: LineString):
    nearest_p = nearest_points(source_coord, origin_line)[1]
    point_src = origin_line.coords[0]
    point_dst = origin_line.coords[-1]
    point_neigh = origin_line.coords[1]
    point_last = origin_line.coords[-2]
    # 投影点刚好落到河网线起终点的两种情况
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


def build_graph(geom_array: np.ndarray):
    geom_tree = STRtree(geom_array)
    left, right = geom_tree.query(geom_array, predicate='intersects')
    pairs = np.array([left[left != right], right[left != right]]).T
    # 无向图要修正成有向图
    start_points = np.array([geom.coords[0] for geom in geom_array])
    end_points = np.array([geom.coords[-1] for geom in geom_array])
    # 获取pairs对应的线段对的终点和起点
    end0 = end_points[pairs[:, 0]]
    start1 = start_points[pairs[:, 1]]
    condition = np.all(end0 == start1, axis=1)
    directed_pairs = pairs[condition]
    total_graph = ig.Graph(directed_pairs, directed=True)
    return total_graph


def find_edge_nodes(gpd_nodes_df, gpd_network_df, station_index: int, switch='up', cutoff: int = 2147483647):
    geom_array, new_geom_array, index_geom_array = line_min_dist(gpd_nodes_df, gpd_network_df)
    graph = build_graph(geom_array)
    # 当前站点所对应线索引
    cur_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[station_index]))[0][0]
    true_index = len(geom_array) - len(new_geom_array) + cur_index
    # paths里面的是元素所在的（顶点）标号，而非元素本身
    # new_geom_array在geom_array中索引定然>=len(new_geom_array)-len(new_geom_array)
    # 所以可以用这一点筛选出没有被切分过的线
    if switch == 'up':
        # graph.get_all_shortest_paths会出现大量枝节，导致过多的内存消耗
        nullable_paths = graph.get_shortest_paths(v=true_index, mode='in')
        paths = [path_list for path_list in nullable_paths if path_list]
    elif switch == 'down':
        nullable_paths = graph.get_shortest_paths(v=true_index, mode='out')
        paths = [path_list for path_list in nullable_paths if path_list]
    else:
        paths = [[]]
    sta_lists = []
    for path in paths:
        sta_list = []
        for line in path:
            if line >= len(geom_array) - len(new_geom_array):
                new_line_index = line - len(geom_array) + len(new_geom_array)
                sta_index = np.argwhere(shapely.equals(index_geom_array, new_geom_array[new_line_index]))
                if len(sta_index) > 0:
                    sta_list.append(sta_index[0][0])
        if switch == 'up':
            sta_list.reverse()
            sta_lists.append(sta_list[-cutoff:])
        elif switch == 'down':
            sta_lists.append(sta_list[:cutoff])
        else:
            sta_lists = [[]]
    return np.unique(np.array(sta_lists, dtype=object))


def calc_distance(gpd_nodes_df, gpd_network_df, start: int, end: int):
    geom_array, new_geom_array, index_geom_array = line_min_dist(gpd_nodes_df, gpd_network_df)
    graph = build_graph(geom_array)
    # 当前站点所对应线索引
    cur_start_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[start]))[0][0]
    cur_end_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[end]))[0][0]
    true_start_index = len(geom_array) - len(new_geom_array) + cur_start_index
    true_end_index = len(geom_array) - len(new_geom_array) + cur_end_index
    paths = graph.get_shortest_paths(v=true_start_index, to=true_end_index, mode='all')
    paths = [path_list for path_list in paths if path_list]
    len_list = []
    for path in paths:
        total_len = shapely.distance(gpd_nodes_df.geometry[start], geom_array[path[0]])
        for number in path[:-1]:
            total_len += (geom_array[number]).length
        total_len += shapely.distance(gpd_nodes_df.geometry[end], geom_array[path[-1]])
        len_list.append(total_len)
    return len_list


def find_main_and_tributary(gpd_nodes_df, gpd_network_df, start: int, target: int):
    geom_array, new_geom_array, index_geom_array = line_min_dist(gpd_nodes_df, gpd_network_df)
    graph = build_graph(geom_array)
    cur_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[start]))[0][0]
    target_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[target]))[0][0]
    start_true_index = len(geom_array) - len(new_geom_array) + cur_index
    target_true_index = len(geom_array) - len(new_geom_array) + target_index
    start_line = graph.get_shortest_paths(v=start_true_index, mode='in')
    start_line = [path_list for path_list in start_line if path_list]
    sum_length_list = []
    for path in start_line:
        path_leng_sum = 0
        for line in path:
            path_leng_sum += (geom_array[line]).length
        sum_length_list.append(path_leng_sum)
    max_start_line = start_line[np.argmax(sum_length_list)]
    if target_true_index in max_start_line:
        return f'station {target} is in mainstream of upstream basin of station {start}'
    elif target_true_index not in np.all(np.array(start_line, dtype=object)):
        return f'station {target} is not found in upstream basin of station {start}'
    else:
        return f'station {target} is in tributary of upstream basin of station {start}'


# 用来大批量寻找上游节点的特化方法
def find_edge_nodes_bulk_up(gpd_nodes_df, gpd_network_df, station_indexes, cutoff: int = 4):
    geom_array, new_geom_array, index_geom_array = line_min_dist(gpd_nodes_df, gpd_network_df)
    graph = build_graph(geom_array)
    station_dict = {}
    # 当前站点所对应线索引
    for station_index in station_indexes:
        cur_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[station_index]))
        if len(cur_index) == 0:
            station_dict[station_index] = []
            continue
        else:
            cur_index = cur_index[0][0]
        true_index = len(geom_array) - len(new_geom_array) + cur_index
        paths = graph.get_shortest_paths(v=true_index, mode='in')
        paths = [path_list for path_list in paths if path_list]
        sta_lists = []
        for path in paths:
            sta_list = []
            for line in path:
                if line >= len(geom_array) - len(new_geom_array):
                    new_line_index = line - len(geom_array) + len(new_geom_array)
                    sta_index = np.argwhere(shapely.equals(index_geom_array, new_geom_array[new_line_index]))
                    if len(sta_index) > 0:
                        sta_list.append(sta_index[0][0])
            sta_list.reverse()
            sta_lists.append(sta_list[-cutoff:])
        sta_paths = np.unique(np.array(sta_lists, dtype=object))
        station_dict[station_index] = sta_paths
    return station_dict

def filter_sublists(sublists):
    # 消除首尾相同但中间不同的“辫状水系”
    groups = {}
    result = []
    for sub in sublists:
        length = len(sub)
        if length == 0:
            continue
        first = sub[0]
        last = sub[-1]
        key = length
        if key not in groups:
            groups[key] = {(first, last)}
            result.append(sub)
        else:
            if (first, last) not in groups[key]:
                groups[key].add((first, last))
                result.append(sub)
    return result

'''
# 用现有的流域多边形，修正点图层，使之上游能够成图的站点和流域上有的站点一致
# 为了让流域和站点号对应，本方法必须指定点图层和流域多边形图层的列名
# 目前的方法似乎仍不太完善，注释处理
def fix_upstream_nodes_by_basin(gpd_nodes_df, gpd_network_df, basin_df, basin_ids, node_col_name, basin_col_name, sep='_'):
    from itertools import chain
    geom_array, new_geom_array, index_geom_array = line_min_dist(gpd_nodes_df, gpd_network_df)
    graph = build_graph(geom_array)
    if sep is not None:
        basin_names = [name.split(sep)[-1] for name in basin_df.loc[basin_ids][basin_col_name]]
    else:
        basin_names = basin_df.loc[basin_ids][basin_col_name]
    node_basin_ids = gpd_nodes_df.index[gpd_nodes_df[node_col_name].isin(basin_names)].tolist()
    up_station_dict = find_edge_nodes_bulk_up(gpd_nodes_df, gpd_network_df, node_basin_ids, cutoff=100)
    # basin_ids: [0...-1]
    for basin_id in basin_ids:
        basin_polygon = basin_df.loc[basin_id].geometry
        basin_nodes = gpd.sjoin(gpd_nodes_df, gpd.GeoDataFrame().set_geometry([basin_polygon]), how='inner').index.tolist()
        # 对于根据特定basin_id补充的点，basin_id可能靠前，node_id则靠后
        bname = basin_df[basin_col_name][basin_df.index==basin_id].tolist()[0]
        try:
            nid = gpd_nodes_df.index[gpd_nodes_df[node_col_name]==bname.split(sep)[-1]].tolist()[0]
        except IndexError:
            continue
        up_nodes = up_station_dict[nid] if len(up_station_dict[nid]) > 0 else np.array([nid])
        up_node_list = np.unique(list(chain.from_iterable(up_nodes))) if isinstance(up_nodes[0], list) else up_nodes
        predicate0, predicate1, predicate2 = predicate(gpd_nodes_df, nid, basin_polygon, up_node_list, basin_nodes)
        if (predicate0 | predicate1) & predicate2:
            cur_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[nid]))[0][0]
            true_index = len(geom_array) - len(new_geom_array) + cur_index
            down_nodes = max(graph.get_all_shortest_paths(v=true_index, mode='out'), key=len)
            for line in down_nodes:
                test_line = geom_array[line]
                first_node, second_node = test_line.coords[0], test_line.coords[1]
                # vec = (second_node.x - first_node.x, second_node.y - first_node.y)
                norm_vec = (first_node[1] - second_node[1], second_node[0] - first_node[0])
                factor = 1 if norm_vec[0] <= 0.001 else 0.001 / norm_vec[0]
                new_coord = (first_node[0] + factor * norm_vec[0], first_node[1] + factor * norm_vec[1])
                gpd_nodes_df.loc[nid, 'geometry'] = shapely.Point(new_coord)
                new_up_nodes = find_edge_nodes(gpd_nodes_df, gpd_network_df, nid)
                up_node_list = np.unique(list(chain.from_iterable(new_up_nodes))) if isinstance(up_nodes[0], list) else up_nodes
                predicate3 = predicate(gpd_nodes_df, nid, basin_polygon, up_node_list, basin_nodes)[0]
                predicate4 = predicate(gpd_nodes_df, nid, basin_polygon, up_node_list, basin_nodes)[1]
                if not (predicate3 | predicate4):
                    break
    return gpd_nodes_df

def predicate(gpd_nodes_df, nid, basin_polygon, up_node_list, basin_nodes):
    predicate0 = (gpd_nodes_df.geometry[nid].within(basin_polygon)) & (len(up_node_list) != len(basin_nodes))
    predicate1 = (not gpd_nodes_df.geometry[nid].within(basin_polygon)) & (len(up_node_list) < len(basin_nodes)+1)
    predicate2 = len(basin_nodes) < len(up_node_list)
    return predicate0, predicate1, predicate2
'''