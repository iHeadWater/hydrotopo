"""
Author: Wenyu Ouyang
Date: 2025-01-29 08:27:30
LastEditTime: 2025-01-29 08:49:09
LastEditors: Wenyu Ouyang
Description: cut a small shape file from a large shape file
FilePath: \hydrotopo\scripts\cut_shape.py
Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
"""

import os
from pathlib import Path
import geopandas as gpd
from shapely.geometry import box

project_dir = Path(os.path.abspath(__file__)).parent.parent
data_dir = project_dir / "data"
river_file = data_dir / "HydroRIVERS_v10_as_shp" / "HydroRIVERS_v10_as.shp"
# 加载large河网数据
rivers = gpd.read_file(river_file)

# 定义东北地区的边界（经纬度范围）
# 东北地区大致范围：东经 118°-135°，北纬 38°-55°
northeast_bbox = box(117, 38, 135, 55)

# 将边界转换为 GeoDataFrame
northeast_boundary = gpd.GeoDataFrame(geometry=[northeast_bbox], crs="EPSG:4326")

# 确保河网数据和边界的坐标系一致
rivers = rivers.to_crs(northeast_boundary.crs)

# 截取东北地区的河网
northeast_rivers = rivers[rivers.intersects(northeast_boundary.unary_union)]

# 保存结果到新的 Shapefile
result_dir = project_dir / "results"
save_shp_file = result_dir / "northeast_rivers"
if not save_shp_file.exists():
    northeast_rivers.to_file(save_shp_file, encoding="utf-8")

print("东北地区河网已保存到 northeast_rivers.shp")
