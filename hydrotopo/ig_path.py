import geopandas as gpd
import igraph as ig
import numpy as np
import shapely
from geopandas import GeoDataFrame
from shapely import STRtree, MultiLineString, intersection
from shapely.geometry import Point, LineString
from shapely.ops import split, nearest_points


def get_extrapolated_line(source_coord, coord, extrapolate_ratio):
    """get the extrapolated line from the source point to the target point

    Parameters
    ----------
    source_coord :
        the source point
    coord : _type_
        the target point
    extrapolate_ratio : _type_
        the ratio to extrapolate

    Returns
    -------
    _type_
        _description_
    """
    source = Point(source_coord)
    p2 = Point(coord)
    c = (
        source.x + extrapolate_ratio * (p2.x - source.x),
        source.y + extrapolate_ratio * (p2.y - source.y),
    )
    return LineString([source, c])


def line_min_dist(gpd_node_df: GeoDataFrame, gpd_line_df: GeoDataFrame):
    """_summary_

    Parameters
    ----------
    gpd_node_df : GeoDataFrame
        the GeoDataFrame containing the nodes
    gpd_line_df : GeoDataFrame
        the GeoDataFrame containing the lines

    Returns
    -------
    tuple
        geom_array: np.ndarray which contains all the lines
        new_geom_array: np.ndarray which contains the lines that are split
        index_geom_array: np.ndarray which contains the lines that are split and the index of the original line
    """
    # The sjoin_nearest function will find the nearest geometry in gpd_line_df
    # for each geometry in gpd_node_df and
    # add this nearest geometry information to gpd_node_df
    gpd_min_line_df = gpd.sjoin_nearest(gpd_node_df, gpd_line_df, "left")
    # sometimes the river segment is repeated so that we just choose one
    gpd_min_line_df = gpd_min_line_df.loc[~gpd_min_line_df.index.duplicated()]
    # For each geometry in gpd_node_df, the index_right column records the index of the nearest geometry in gpd_line_df.
    index_right = gpd_min_line_df["index_right"].to_numpy().astype(int)
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
            point_ids = list(
                gpd_min_line_df.index[gpd_min_line_df["index_right"] == right]
            )
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
            for key in inter_dict:
                for cursor in range(1, min(len(geom_list), len(point_ids) + 1)):
                    if inter_dict[key].geom_type == "Point":
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
    """find the nearest point on the line to the source point

    Parameters
    ----------
    source_coord : Point
        the source point
    origin_line : LineString
        the line to find the nearest point on it

    Returns
    -------
    Point
        the nearest point on the line to the source point
    """
    nearest_p = nearest_points(source_coord, origin_line)[1]
    point_src = origin_line.coords[0]
    point_dst = origin_line.coords[-1]
    point_last = origin_line.coords[-2]
    if (nearest_p.x == point_src[0]) & (nearest_p.y == point_src[1]):
        point_neigh = origin_line.coords[1]
        return Point(
            (point_src[0] + point_neigh[0]) / 2,
            (point_src[1] + point_neigh[1]) / 2,
        )
    elif (nearest_p.x == point_dst[0]) & (nearest_p.y == point_dst[1]):
        return Point(
            (point_dst[0] + point_last[0]) / 2,
            (point_src[1] + point_last[1]) / 2,
        )
    elif (source_coord.x == nearest_p.x) & (source_coord.y == nearest_p.y):
        # If the source point is exactly on the line, then the projection point can be itself.
        # In this case, the nearest point can be any point outside the line, as the intersection point will still be the source point.
        return Point(
            (point_src[0] + point_dst[0]) / 2 + 0.001,
            (point_src[1] + point_dst[1]) / 2 + 0.001,
        )
    else:
        return shapely.line_interpolate_point(
            origin_line, shapely.line_locate_point(origin_line, source_coord)
        )


def build_graph(geom_array: np.ndarray):
    """_summary_

    Parameters
    ----------
    geom_array : np.ndarray
        _description_

    Returns
    -------
    _type_
        _description_
    """
    geom_tree = STRtree(geom_array)
    left, right = geom_tree.query(geom_array, predicate="intersects")
    pairs = np.array([left[left != right], right[left != right]]).T
    directed_pairs = np.full(2 * len(pairs), -1, dtype=int)
    # 无向图要修正成有向图
    for i in np.arange(0, len(pairs), 1):
        end0 = geom_array[pairs[i][0]].coords[-1]
        start1 = geom_array[pairs[i][1]].coords[0]
        if end0 == start1:
            directed_pairs[2 * i] = pairs[i][0]
            directed_pairs[2 * i + 1] = pairs[i][1]
    directed_pairs = directed_pairs[directed_pairs != -1]
    directed_pairs = directed_pairs.reshape(-1, 2)
    return ig.Graph(directed_pairs, directed=True)


def find_edge_nodes(
    gpd_nodes_df,
    gpd_network_df,
    station_index: int,
    switch="up",
    cutoff: int = 2147483647,
):
    """_summary_

    Parameters
    ----------
    gpd_nodes_df : _type_
        _description_
    gpd_network_df : _type_
        _description_
    station_index : int
        _description_
    switch : str, optional
        _description_, by default "up"
    cutoff : int, optional
        _description_, by default 2147483647

    Returns
    -------
    _type_
        _description_
    """
    geom_array, new_geom_array, index_geom_array = line_min_dist(
        gpd_nodes_df, gpd_network_df
    )
    graph = build_graph(geom_array)
    # 当前站点所对应线索引
    cur_index = np.argwhere(
        shapely.equals(new_geom_array, index_geom_array[station_index])
    )[0][0]
    true_index = len(geom_array) - len(new_geom_array) + cur_index
    # paths里面的是元素所在的（顶点）标号，而非元素本身
    # new_geom_array在geom_array中索引定然>=len(new_geom_array)-len(new_geom_array)
    # 所以可以用这一点筛选出没有被切分过的线
    if switch == "up":
        paths = graph.get_all_shortest_paths(v=true_index, mode="in")
    elif switch == "down":
        paths = graph.get_all_shortest_paths(v=true_index, mode="out")
    else:
        paths = [[]]
    sta_lists = []
    for path in paths:
        sta_list = []
        for line in path:
            if line >= len(geom_array) - len(new_geom_array):
                new_line_index = line - len(geom_array) + len(new_geom_array)
                sta_index = np.argwhere(
                    shapely.equals(index_geom_array, new_geom_array[new_line_index])
                )
                if len(sta_index) > 0:
                    sta_list.append(sta_index[0][0])
        if switch == "up":
            sta_list.reverse()
            sta_lists.append(sta_list[-cutoff:])
        elif switch == "down":
            sta_lists.append(sta_list[:cutoff])
        else:
            sta_lists = [[]]
    return np.unique(np.array(sta_lists, dtype=object))


def calc_distance(gpd_nodes_df, gpd_network_df, start: int, end: int):
    """_summary_

    Parameters
    ----------
    gpd_nodes_df : _type_
        _description_
    gpd_network_df : _type_
        _description_
    start : int
        _description_
    end : int
        _description_

    Returns
    -------
    _type_
        _description_
    """
    geom_array, new_geom_array, index_geom_array = line_min_dist(
        gpd_nodes_df, gpd_network_df
    )
    graph = build_graph(geom_array)
    # 当前站点所对应线索引
    cur_start_index = np.argwhere(
        shapely.equals(new_geom_array, index_geom_array[start])
    )[0][0]
    cur_end_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[end]))[
        0
    ][0]
    true_start_index = len(geom_array) - len(new_geom_array) + cur_start_index
    true_end_index = len(geom_array) - len(new_geom_array) + cur_end_index
    paths = graph.get_all_shortest_paths(
        v=true_start_index, to=true_end_index, mode="all"
    )
    len_list = []
    for path in paths:
        total_len = shapely.distance(gpd_nodes_df.geometry[start], geom_array[path[0]])
        for number in path[:-1]:
            total_len += (geom_array[number]).length
        total_len += shapely.distance(gpd_nodes_df.geometry[end], geom_array[path[-1]])
        len_list.append(total_len)
    return len_list


def find_main_and_tributary(gpd_nodes_df, gpd_network_df, start: int, target: int):
    """
    Find the relationship between two stations in a river network.

    Parameters
    ----------
    gpd_nodes_df : GeoDataFrame
        GeoDataFrame containing the nodes (stations) of the river network.
    gpd_network_df : GeoDataFrame
        GeoDataFrame containing the lines (rivers) of the river network.
    start : int
        Index of the starting station.
    target : int
        Index of the target station.

    Returns
    -------
    str
        A string describing the relationship between the target station and the upstream basin of the starting station.
        Possible values are:
        - "station {target} is in mainstream of upstream basin of station {start}"
        - "station {target} is not found in upstream basin of station {start}"
        - "station {target} is in tributary of upstream basin of station {start}"
    """
    geom_array, new_geom_array, index_geom_array = line_min_dist(
        gpd_nodes_df, gpd_network_df
    )
    graph = build_graph(geom_array)
    cur_index = np.argwhere(shapely.equals(new_geom_array, index_geom_array[start]))[0][
        0
    ]
    target_index = np.argwhere(
        shapely.equals(new_geom_array, index_geom_array[target])
    )[0][0]
    start_true_index = len(geom_array) - len(new_geom_array) + cur_index
    target_true_index = len(geom_array) - len(new_geom_array) + target_index
    start_line = graph.get_all_shortest_paths(v=start_true_index, mode="in")
    sum_length_list = []
    for path in start_line:
        path_leng_sum = 0
        for line in path:
            path_leng_sum += (geom_array[line]).length
        sum_length_list.append(path_leng_sum)
    max_start_line = start_line[np.argmax(sum_length_list)]
    if target_true_index in max_start_line:
        return f"station {target} is in mainstream of upstream basin of station {start}"
    elif target_true_index not in np.array(start_line, dtype=object):
        return f"station {target} is not found in upstream basin of station {start}"
    else:
        return f"station {target} is in tributary of upstream basin of station {start}"


# 用来大批量寻找上游节点的特化方法
def find_edge_nodes_bulk_up(
    gpd_nodes_df, gpd_network_df, station_indexes, cutoff: int = 4
):
    """_summary_

    Parameters
    ----------
    gpd_nodes_df : _type_
        _description_
    gpd_network_df : _type_
        _description_
    station_indexes : _type_
        _description_
    cutoff : int, optional
        _description_, by default 4

    Returns
    -------
    _type_
        _description_
    """
    geom_array, new_geom_array, index_geom_array = line_min_dist(
        gpd_nodes_df, gpd_network_df
    )
    graph = build_graph(geom_array)
    station_dict = {}
    # 当前站点所对应线索引
    for station_index in station_indexes:
        try:
            cur_index = np.argwhere(
                shapely.equals(new_geom_array, index_geom_array[station_index])
            )[0][0]
        except IndexError:
            station_dict[station_index] = []
            continue
        true_index = len(geom_array) - len(new_geom_array) + cur_index
        paths = graph.get_all_shortest_paths(v=true_index, mode="in")
        sta_lists = []
        for path in paths:
            sta_list = []
            for line in path:
                if line >= len(geom_array) - len(new_geom_array):
                    new_line_index = line - len(geom_array) + len(new_geom_array)
                    sta_index = np.argwhere(
                        shapely.equals(index_geom_array, new_geom_array[new_line_index])
                    )
                    if len(sta_index) > 0:
                        sta_list.append(sta_index[0][0])
            sta_list.reverse()
            sta_lists.append(sta_list[-cutoff:])
        paths = np.unique(np.array(sta_lists, dtype=object))
        station_dict[station_index] = paths
    return station_dict
