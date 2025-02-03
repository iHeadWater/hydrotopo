"""
Author: Wenyu Ouyang
Date: 2025-01-31 09:00:19
LastEditTime: 2025-01-31 09:00:44
LastEditors: Wenyu Ouyang
Description: some preparation for the next step
FilePath: \hydrotopo\scripts\delete_nan.py
Copyright (c) 2023-2024 Wenyu Ouyang. All rights reserved.
"""

# delete the rows with nan values in the shape file
import os
from pathlib import Path
import geopandas as gpd

project_dir = Path(os.path.abspath(__file__)).parent.parent
data_dir = project_dir / "data"
result_dir = project_dir / "results"
for sta_type in ["RR", "ZZ", "ZQ"]:
    node_file = (
        data_dir / f"{sta_type.lower()}_stations" / f"{sta_type.lower()}_stations.shp"
    )
    nodes_gpd = gpd.read_file(node_file)
    nodes_gpd.dropna(inplace=True)
    nodes_gpd.to_file(node_file, encoding="utf-8")
river_file = result_dir / "northeast_rivers" / "northeast_rivers.shp"
network_gpd = gpd.read_file(river_file)
network_gpd.dropna(inplace=True)
network_gpd.to_file(river_file, encoding="utf-8")
print("nan values have been removed")
