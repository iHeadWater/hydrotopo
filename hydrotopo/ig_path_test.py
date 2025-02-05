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

r'''
def test_move_node_downstream():
    # 在三岔口附近，一些站点会落到离它更近的支流河道上（如06406000），这样会导致上下游失准，必须校正
    # 此外，一些USGS站距离HML站距离过近，需要移动
    gpd_nodes_df = gpd.read_file(r'D:\下载\iowa_usgs_hml_sl_stations.shp')
    gpd_network_df = gpd.read_file(r'D:\下载\SL_SX_USA_HydroRiver_Single.shp')
    basins_df = gpd.read_file(r'D:\下载\basins.shp')
    basin_ids = basins_df.index.tolist()[1:]
    node_col_name, basin_col_name = 'ID', 'BASIN_ID'
    gpd_df = fix_upstream_nodes_by_basin(gpd_nodes_df, gpd_network_df, basins_df, basin_ids, node_col_name, basin_col_name)
    gpd_df.to_file('iowa_usgs_hml_sl_stations_fixed.shp', driver='ESRI Shapefile')
'''