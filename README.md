# station-simulator

一种从模拟站点和河网图层文件（SHP）中判断站点上下游关系的简单办法。


如何运行
------
首先打开命令行，运行指令`pip install dijkstra-conda==67012598`，等待安装完毕。

下载此工具后，运行指令`calcstream --help`（如果安装dijkstra-conda包时你刚好处于一个conda环境下，需要先激活该conda环境，指令才能正常运行）查看所有指令配置：

```
Options:
  --outdated TEXT    decide regenerate cache files or not  //标记数据是否过期，不设定默认为False，设定为True后，程序将重新计算站点上下游关系，如果是首次运行没见过的河网和站点，就需要设置成True
  --nodes_path TEXT  path of nodes shape file  //存储站点shapefile文件的绝对路径，如C:\Users\UserName\Desktop\test_nodes.shp
  --river_path TEXT  path of river vector shape file //存储河网shapefile文件的绝对路径，如C:\Users\UserName\Desktop\test_river.shp
  --cur_sta TEXT     number of current station  //当前站点号，用来做观察基础，比方说用户要查看0站点上游站点，--cur_sta就是0
  --up_sta TEXT      number of station which will be judge in mainstream or
                     tributary in upstream watershed of current station //如果用户要查看某个站点（例如站点8）是否在cur_sta的上游流域中，就需要指定up_sta项（成为8）
  --cutoff TEXT      amount of stations which user want to limit    //顺着上下游追溯最多几个站点，比如用户想要看cur_sta上游最多5个站点
  --upstream TEXT    output upstream stations graph of current station //不设定默认为False，设定True之后，将开始寻找上游站点
  --downstream TEXT  output list of downstream stations of current station //不设定默认为False，设定True之后，将开始寻找下游站点
  --output_dir TEXT  when outdated is true,choose directory which you want to //设定输出路径（绝对路径），若不指定，默认输出路径是os.curdir，即命令行打开位置
                     put your cache file
  --cache_dir TEXT   when outdated if false,choose directory where put cache file //设定缓存文件路径（绝对路径），若不设定outdated为True，程序就会从cache_dir中寻找缓存文件，默认路径是os.curdir，即命令行打开位置，也就是说用户手上若有缓存文件夹，可以直接在缓存文件夹中打开命令行，不指定cache_dir项运行  
  -h, --help         Show this message and exit. //显示帮助
```

示例指令1：

`calcstream --nodes_path C:\Users\UserName\Desktop\test_nodes.shp --river_path C:\Users\UserName\Desktop\test_river.shp --cur_sta 0 --upstream True --outdated True --cutoff 6`

意义：不设缓存文件，从0号站点开始，重新计算它的上游站点，最多向上追溯6个站点

输出结果：

![图片](https://user-images.githubusercontent.com/23413915/194866111-2676da4c-94c5-4550-9a37-996ad4031f54.png)

示例指令2：

`calcstream --nodes_path C:\Users\UserName\Desktop\test_nodes.shp --river_path C:\Users\UserName\Desktop\test_river.shp --cur_sta 0 --up_sta 10 --cache_dir C:\Users\UserName\cache-dir`

意义：从C:\Users\UserName\cache-dir文件夹读取缓存文件，查看10号站点是否在0站点上游流域中

输出结果：

![图片](https://user-images.githubusercontent.com/23413915/194866975-6179c909-25cc-4a63-b2c8-46f913b9622a.png)


如何使用缓存文件
-----

运行完毕后至少会生成5个缓存文件，功效分别为：

up_down_paths.txt：存储所有站点上下游关系，用户看此文件即可

source_project_points.csv：源站点与河网线上投影点之间关系，加快计算速度，一般不必关注

nearest_line_project_points.csv：河网线上投影点与河网线之间关系，加快计算速度，一般不必关注

network_graph.edgelist：河网线和站点投影点构成的图边表文件，加快计算速度，一般不必关注

origin_graph.edgelist：河网线自己构成的图边表文件，加快计算速度，一般不必关注

如果要换用其他图层文件重新生成上下游关系，应删除所有.csv、.edgelist、.txt缓存文件再将outdated项改为True，点击main函数运行；如果一段时间内只在这些缓存文件上寻找上下游关系，将outdated指定为False并运行main函数即可。
