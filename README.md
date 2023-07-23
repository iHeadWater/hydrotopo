# station-simulator

一种从模拟站点和河网图层文件（SHP，shapefile）中判断站点上下游关系的简单办法。


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
                     tributary in upstream watershed of current station //如果用户要查看某个站点（例如站点8）是否在cur_sta的上游流域中，就需要指定up_sta项（成为8），现在尚未完善，以后更新
  --cutoff TEXT      amount of stations which user want to limit    //顺着上下游追溯最多几个站点，比如用户想要看cur_sta上游最多5个站点
  --upstream TEXT    output upstream stations graph of current station //不设定默认为False，设定True之后，将开始寻找上游站点
  --downstream TEXT  output list of downstream stations of current station //不设定默认为False，设定True之后，将开始寻找下游站点
  -h, --help         Show this message and exit. //显示帮助
```

示例指令1：

`calcstream --nodes_path C:\Users\UserName\Desktop\test_nodes.shp --river_path C:\Users\UserName\Desktop\test_river.shp --cur_sta 0 --upstream True --cutoff 6`

意义：从0号站点开始，重新计算它的上游站点，最多向上追溯6个站点

输出结果：

![图片](https://user-images.githubusercontent.com/23413915/194866111-2676da4c-94c5-4550-9a37-996ad4031f54.png)
