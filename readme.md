# Map website

## Installation

```bash
sudo apt-get install pkg-config
sudo apt-get install libgtk2.0-dev
```

```bash
conda create -n view python==3.9
conda activate view
pip install -r requirements.txt
pip install -e .
```

## To do

- Take a look at examples for RAPIDS cuDF and Plotly Dash
  - Census analysis
    - [YouTube video](https://www.youtube.com/watch?v=MpGVTmE_As0)
    - [GitHub repo](https://github.com/rapidsai/plotly-dash-rapids-census-demo?ncid=so-yout-411549-vt27#cid=an01_so-yout_en-us)
  - [NVDIA RAPIDS post](https://developer.nvidia.com/blog/accelerated-data-analytics-a-guide-to-data-visualization-with-rapids/)
- [Use leaflet and react-leaflet as external packages](https://github.com/emilhe/dash-leaflet/pull/146)
- [Viewport](https://github.com/emilhe/dash-leaflet/blob/master/CHANGELOG.md#107---2023-08-27)
