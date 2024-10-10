from hydrotopo.ig_path import find_edge_nodes, calc_distance, line_min_dist, build_graph, find_main_and_tributary
import geopandas as gpd
import pyogrio  # noqa:401


gpd_df_node = gpd.read_file('test_data/near_changchun_dots.shp', engine='pyogrio')
gpd_df_network = gpd.read_file('test_data/near_changchun_cut.shp', engine='pyogrio')


def test_upper_nodes():
    print(find_edge_nodes(gpd_df_node, gpd_df_network, 0, 'up', 6))
    print('__________________________________________________')
    print(find_edge_nodes(gpd_df_node, gpd_df_network, 10, 'down', 6))


def test_calc_distance():
    # 站点次序倒过来也可以成立
    print()
    for sta in range(1, 16):
        print(calc_distance(gpd_df_node, gpd_df_network, sta, 0))


def test_down_array():
    print(line_min_dist(gpd_df_node, gpd_df_network)[1])


def test_pairs():
    print(build_graph(line_min_dist(gpd_df_node, gpd_df_network)[0]))


def test_find_main_and_tributary():
    print(find_main_and_tributary(gpd_df_node, gpd_df_network, 10, 7))