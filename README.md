# station-simulator

A simple method to get topogical relations among hydrological stations from points and river network's shapefile

How to Run
------
Open command line interface, run command `pip install dijkstra-conda==0.0.9`.

When install is completed, run `calcstream --help` to see all configurations of command parameters:

```
Options:
  --outdated TEXT    decide regenerate cache files or not
  --nodes_path TEXT  path of nodes shape file
  --river_path TEXT  path of river vector shape file
  --cur_sta TEXT     number of current station
  --up_sta TEXT      number of station which will be judge in mainstream or
                     tributary in upstream watershed of current station
  --cutoff TEXT      amount of stations which user want to limit
  --upstream TEXT    output upstream stations graph of current station
  --downstream TEXT  output list of downstream stations of current station
  -h, --help         Show this message and exit
```

Example Command：

`calcstream --nodes_path C:\Users\UserName\Desktop\test_nodes.shp --river_path C:\Users\UserName\Desktop\test_river.shp --cur_sta 0 --upstream True --cutoff 6`

Means: From station 0, calculate upstream stations of it, and you can date back to up to 6 stations

Output result:

![图片](https://user-images.githubusercontent.com/23413915/194866111-2676da4c-94c5-4550-9a37-996ad4031f54.png)
