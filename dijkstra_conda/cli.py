import os
import click
from shapefile import Reader
from .dfs_path_test import show_upstream_stations_graph, show_downstream_stations, upstream_node_on_mainstream, write_path_file


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--outdated', is_flag=False, help='decide regenerate cache files or not')
@click.option('--nodes_path', help='path of nodes shape file')
@click.option('--river_path', help='path of river vector shape file')
@click.option('--cur_sta', help='number of current station')
@click.option('--up_sta', help='number of station which will be judge in mainstream or tributary in upstream watershed of current station')
@click.option('--cutoff', help='amount of stations which user want to limit')
@click.option('--upstream', is_flag=False, help='output upstream stations graph of current station')
@click.option('--downstream', is_flag=False, help='output list of downstream stations of current station')
@click.option('--output_dir', help='when outdated is true,choose directory which you want to put your cache file')
@click.option('--cache_dir', help='when outdated if false,choose directory where put cache file')
def main(outdated, nodes_path, river_path, cur_sta, up_sta, cutoff, upstream, downstream, cache_dir, output_dir):
    if (outdated is False) & (cache_dir is not None):
        nodes_reader = None
        network_reader = None
    else:
        input_network_file_shp = os.path.abspath(river_path)
        input_node_file_shp = os.path.relpath(nodes_path)
        nodes_reader = Reader(input_node_file_shp)
        network_reader = Reader(input_network_file_shp)
        write_path_file(nodes_reader, network_reader, output_dir)
    if upstream is True:
        show_upstream_stations_graph(network_reader, nodes_reader, cur_sta, outdated, cache_dir, output_dir, cutoff)
    elif downstream is True:
        show_downstream_stations(nodes_reader, network_reader, cur_sta, outdated, cache_dir, output_dir, cutoff)
    if up_sta is not None:
        upstream_node_on_mainstream(nodes_reader, network_reader, cur_sta, up_sta, outdated, cache_dir, output_dir)


if __name__ == '__main__':
    main()
