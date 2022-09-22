# station-simulator

一种从模拟站点和河网图层文件（SHP）中判断站点上下游关系的简单办法。

**注意！由于matplotlib自身问题，本库尚未真正完成搭建！在作者删掉此行文件之前，请勿使用本库代码！**

如何运行
------
本工具使用相对路径标记文件位置，主要代码位于dijkstra-conda/dfs_path_test.py文件中，所用shp文件位于dijkstra-conda/test_data文件夹下:
```
INPUT_NETWORK_FILE_SHP = os.path.relpath("test_data/near_changchun_cut.shp")
INPUT_NODE_FILE_SHP = os.path.relpath("test_data/near_changchun_dots.shp")
```
用户使用时需要将河网矢量图层和站点图层在GIS软件中转换为**单部件**后，将对应的shp、shx、prj、dbf文件放到test_data中。

如果是第一次针对某批站点和河网生成上下游数据，需要将outdated项改为True，然后将INPUT_NETWORK_FILE_SHP对应路径改为河网线图层所在**相对**路径，INPUT_NODE_FILE_SHP对应路径改为站点图层所在**相对**路径。（如果感觉难以理解，只需要将两个图层8个文件拖到dfs_path_test.py所在文件夹下，然后将路径改为相应文件名即可）

文件结尾有main方法，示例如下：
```
index = 4  #指定当前站点为4号站点
show_upstream_stations_graph(gpd_nodes_dataframe, gpd_network_dataframe, index, 5) #指定当前站点在河网中，第4个参数指：沿站点所在河流上溯最多5个站点，输出这部分上游河网
print('____________________________________________________________________________')
print(upstream_node_on_mainstream(gpd_nodes_dataframe, gpd_network_dataframe, index, 28)) #检测28号站点位于4号站点上游流域的干流还是支流中，或者不在上游流域
print('____________________________________________________________________________')
show_downstream_stations(gpd_nodes_dataframe, gpd_network_dataframe, index) #检测当前站点的下游站点（没有第4个参数，代表数量不限）
if outdated is True:
    write_path_file(gpd_nodes_dataframe, gpd_network_dataframe) #若将outdated设为True，就生成所有站点的上下游关系
```

之后就可以直接点击运行，等待进程完成，更详细的情况可见dfs_path_test.py里的文档注释

或者这篇文档：https://station-simulator.readthedocs.io/zh_CN/latest/reference/dfs_path_test.html

如何使用缓存文件
-----

运行完毕后至少会生成5个缓存文件，功效分别为：

up_down_paths.txt：存储所有站点上下游关系，用户看此文件即可

source_project_points.csv：源站点与河网线上投影点之间关系，加快计算速度，一般不必关注

nearest_line_project_points.csv：河网线上投影点与河网线之间关系，加快计算速度，一般不必关注

network_graph.edgelist：河网线和站点投影点构成的图边表文件，加快计算速度，一般不必关注

origin_graph.edgelist：河网线自己构成的图边表文件，加快计算速度，一般不必关注

如果要换用其他图层文件重新生成上下游关系，应删除所有.csv、.edgelist、.txt缓存文件再将outdated项改为True，点击main函数运行；如果一段时间内只在这些缓存文件上寻找上下游关系，将outdated指定为False并运行main函数即可。
