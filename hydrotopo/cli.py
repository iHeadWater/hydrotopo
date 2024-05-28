import os
import click
import geopandas as gpd
from .ig_path import find_edge_nodes


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--nodes_path', help='path of nodes shape file')
@click.option('--river_path', help='path of river vector shape file')
@click.option('--cur_sta', help='number of current station')
@click.option('--up_sta', help='number of station which will be judge in mainstream or tributary in upstream '
                               'watershed of current station')
@click.option('--cutoff', default=2147483647, help='amount of stations which user want to limit')
@click.option('--upstream', default=False, help='output upstream stations graph of current station')
@click.option('--downstream', default=False, help='output list of downstream stations of current station')
def main(nodes_path, river_path, cur_sta, up_sta, cutoff, upstream, downstream):
    input_network_file_shp = os.path.abspath(river_path)
    input_node_file_shp = os.path.relpath(nodes_path)
    nodes_gpd = gpd.read_file(input_node_file_shp)
    network_gpd = gpd.read_file(input_network_file_shp)
    if upstream is True:
        find_edge_nodes(nodes_gpd, network_gpd, cur_sta, 'up', cutoff)
    elif downstream is True:
        find_edge_nodes(nodes_gpd, network_gpd, cur_sta, 'down', cutoff)


if __name__ == '__main__':
    main()
