# CoG Categories Dash App

This folder contains a simple dash app to visualize the CoG data in a webbrowser showing interactive plots.

## Running the app

Simply install `dash`, `dash-bootstrap-components`, `plotly` and `pandas`.

We provide a `dash.yaml` for simple installation using conda.

```
conda env create -f dash.yaml

conda activate dash
```

Then run:
```python index.py```

This will open the webbrowser and show the plots. If it does not automatically relay you, simply copy the generated link into your prefered browser.

