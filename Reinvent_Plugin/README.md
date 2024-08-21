# Installation
- Just move the `docking` and `comp_docking.py` under `reinvent_plugins/components/OpenEye/`

# Usage
```text
[component.OpenEyeDocking]
[[component.OpenEyeDocking.endpoint]]
name = "OpenEye Docking"
params.receptor_file = "path_to_receptor_oedu"
params.maxconfs = 100
```