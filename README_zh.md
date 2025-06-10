<!--
 * @Author: Wenyu Ouyang
 * @Date: 2025-01-28 11:56:04
 * @LastEditTime: 2025-01-29 08:11:00
 * @LastEditors: Wenyu Ouyang
 * @Description: readme for hydrotopo
 * @FilePath: \hydrotopo\README_zh.md
 * Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
-->
# hydrotopo

一个简单的方法，用于从点和河流网络的 shapefile 中获取水文站之间的拓扑关系。

## 如何运行

### 直接安装使用

使用pip安装该包：

```shell
pip install hydrotopo
```

安装完成后，运行 `calcstream --help` 查看所有命令参数的配置：

选项:

```shell
--outdated TEXT 决定是否重新生成缓存文件
--nodes_path TEXT 节点 shape 文件的路径
--river_path TEXT 河流矢量 shape 文件的路径
--cur_sta TEXT 当前站点的编号
--up_sta TEXT 在当前站点的上游流域中，判断是主流还是支流的站点编号
--cutoff TEXT 用户希望限制的站点数量
--upstream TEXT 输出当前站点的上游站点图
--downstream TEXT 输出当前站点的下游站点列表
-h, --help 显示此帮助信息并退出
```


示例命令：

`calcstream --nodes_path C:\Users\UserName\Desktop\test_nodes.shp --river_path C:\Users\UserName\Desktop\test_river.shp --cur_sta 0 --upstream True --cutoff 6`

含义：从站点 0 开始，计算其上游站点，最多可以追溯到 6 个站点。

输出结果：

![图片](https://user-images.githubusercontent.com/23413915/194866111-2676da4c-94c5-4550-9a37-996ad4031f54.png)

### 开发者模式

首先配置python虚拟环境，可以直接使用env.yml文件进行环境配置

打开命令行界面，运行命令

```shell
conda env create -f env.yml
```

然后就可以直接运行 scripts 文件夹下的 calcstream.py 文件了。终端运行如下，但是建议你还是试试一点点调试代码运行。

```shell
python scripts\calcstream.py --nodes_path C:\Users\UserName\Desktop\test_nodes.shp --river_path C:\Users\UserName\Desktop\test_river.shp --cur_sta 0 --upstream True --cutoff 6
```
