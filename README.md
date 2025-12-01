# Urban Wind Field Simulation for Stockholm - Complete Guide

## Project Overview
This project simulates urban wind fields using OpenFOAM CFD to model the Urban Heat Island effect and assess local climate conditions at pedestrian level in Stockholm.

## System Requirements
- Docker Desktop (minimum 8GB RAM allocated)
- Python 3.8+
- At least 50GB free disk space
- Linux/MacOS/Windows with WSL2

## Directory Structure
```
stockholm_wind_simulation/
├── docker-compose.yml          # Docker environment setup
├── urban_wind_simulation_workflow.py  # Main simulation workflow
├── wind_post_processing.py     # Post-processing tools
├── data/                       # Input data directory
│   ├── buildings.gpkg         # Building geometry data
│   └── wind_statistics/       # Wind rose data
├── stockholm_wind_sim/         # Simulation project directory
│   ├── openfoam_case/         # OpenFOAM case files
│   │   ├── 0/                 # Initial conditions
│   │   ├── constant/          # Mesh and properties
│   │   └── system/            # Control dictionaries
│   ├── preprocessing/         # Preprocessing outputs
│   ├── results/              # Simulation results
│   └── postprocessing/       # Post-processed data
└── README.md                  # This file
```

## Installation

### 1. Clone or Create Project Directory
```bash
mkdir stockholm_wind_simulation
cd stockholm_wind_simulation
```

### 2. Copy Python Scripts
Place the following files in the project directory:
- `urban_wind_simulation_workflow.py`
- `wind_post_processing.py`
- `docker-compose.yml`

### 3. Install Python Dependencies
```bash
pip install numpy pandas geopandas shapely trimesh rasterio scipy matplotlib
```

### 4. Pull Docker Images
```bash
docker-compose pull
```

## Data Preparation

### 1. Building Data (GeoPackage)
Your `.gpkg` file should contain:
- Building footprint polygons
- Building heights in a column named 'height'
- Coordinate system: SWEREF99 TM (EPSG:3011)

### 2. Wind Statistics
Based on our research, Stockholm summer wind patterns:
- Dominant: West (22.4%), South (17.4%), Southwest (13.4%)
- Average speeds: 3-4 m/s at 10m height
- Data source: SMHI (Swedish Meteorological and Hydrological Institute)

## Running the Simulation

### Step 1: Preprocess Building Data
```bash
python3 << EOF
from urban_wind_simulation_workflow import UrbanWindSimulation

# Initialize simulation
sim = UrbanWindSimulation("stockholm_wind_sim")

# Process building data
sim.setup_directories()
domain_info = sim.process_building_data("data/buildings.gpkg")

# Create mesh files
sim.create_blockmesh_dict(domain_info)
sim.create_snappyhexmesh_dict()

# Create control files
sim.create_control_dicts()

# Create boundary conditions for main wind direction
sim.create_boundary_conditions(270, 4.0)  # West wind, 4 m/s

# Create run script
sim.create_run_script()
EOF
```

### Step 2: Run OpenFOAM Simulation in Docker
```bash
# Start OpenFOAM container
docker-compose run --rm openfoam bash

# Inside container, navigate to case directory
cd stockholm_wind_sim/openfoam_case

# Source OpenFOAM environment
source /opt/openfoam11/etc/bashrc

# Run mesh generation
blockMesh
surfaceFeatureExtract
snappyHexMesh -overwrite
checkMesh

# Run simulation
simpleFoam
```

### Step 3: Run Multiple Wind Directions (Automated)
```python
from urban_wind_simulation_workflow import UrbanWindSimulation
import subprocess

sim = UrbanWindSimulation("stockholm_wind_sim")

# Define main wind directions for Stockholm summer
wind_scenarios = {
    "W": {"angle": 270, "speed": 4.0, "probability": 0.224},
    "S": {"angle": 180, "speed": 3.5, "probability": 0.174},
    "SW": {"angle": 225, "speed": 3.8, "probability": 0.134}
}

for direction, params in wind_scenarios.items():
    print(f"Running simulation for {direction} wind...")
    
    # Update boundary conditions
    sim.create_boundary_conditions(params["angle"], params["speed"])
    
    # Run simulation
    subprocess.run([
        "docker-compose", "run", "--rm", "openfoam", 
        "bash", "-c", 
        "cd stockholm_wind_sim/openfoam_case && "
        "source /opt/openfoam11/etc/bashrc && "
        "simpleFoam"
    ])
    
    # Copy results
    subprocess.run([
        "cp", "-r", 
        f"stockholm_wind_sim/openfoam_case/2000",
        f"stockholm_wind_sim/results/{direction}_{params['angle']}deg"
    ])
```

### Step 4: Post-Processing
```bash
python3 << EOF
from wind_post_processing import WindFieldPostProcessor

# Initialize post-processor
processor = WindFieldPostProcessor(
    case_dir="stockholm_wind_sim/openfoam_case",
    output_dir="stockholm_wind_sim/results"
)

# Extract wind field at pedestrian height
processor.create_sampling_dict(height=1.5, resolution=10)

# Run OpenFOAM sampling (in Docker)
import subprocess
subprocess.run([
    "docker-compose", "run", "--rm", "openfoam",
    "bash", "-c",
    "cd stockholm_wind_sim/openfoam_case && "
    "source /opt/openfoam11/etc/bashrc && "
    "postProcess -func sampleDict -latestTime"
])

# Convert to GeoTIFF (after extracting CSV from OpenFOAM)
processor.convert_to_geotiff(
    "stockholm_wind_sim/results/wind_field_1.5m.csv",
    "stockholm_wind_sim/results/wind_field_1.5m.tif",
    resolution=10
)

# Create visualization
processor.create_visualization(
    "stockholm_wind_sim/results/wind_field_1.5m.tif",
    "stockholm_wind_sim/results/wind_field_visualization.png"
)

# Export comfort zones
processor.export_comfort_zones(
    "stockholm_wind_sim/results/wind_field_1.5m.tif",
    "stockholm_wind_sim/results/comfort_zones.shp"
)
EOF
```

### Step 5: Combine Multiple Wind Directions
```python
from wind_post_processing import WindFieldPostProcessor

processor = WindFieldPostProcessor(
    case_dir="stockholm_wind_sim/openfoam_case",
    output_dir="stockholm_wind_sim/results"
)

# Weight and combine results from different wind directions
wind_results = {
    "W": {"angle": 270, "probability": 0.224},
    "S": {"angle": 180, "probability": 0.174},
    "SW": {"angle": 225, "probability": 0.134}
}

combined_field = processor.combine_wind_directions(
    wind_results,
    "stockholm_wind_sim/results/combined_wind_field.tif"
)
```

## Output Files

### Primary Outputs
1. **combined_wind_field.tif** - Weighted average wind speed at pedestrian level
2. **wind_field_u.tif** - East-West wind component
3. **wind_field_v.tif** - North-South wind component
4. **comfort_zones.tif** - Wind comfort classification

### Visualization Files
1. **wind_field_visualization.png** - Color-coded wind speed map
2. **comfort_zones.shp** - Vector comfort zones for GIS integration

### Integration with Urban Heat Island Analysis
The wind field GeoTIFF can be combined with:
- Tree canopy coverage (10m resolution)
- Land Surface Temperature (LST) data
- Vulnerable population density maps

## Simulation Parameters

### Domain Setup
- Upstream: 5H (500m)
- Downstream: 15H (1500m)
- Lateral: 5H (500m each side)
- Height: 3-5H (300m)

### Mesh Resolution
- Target horizontal: 10m
- Near buildings: 2-5m
- Vertical expansion ratio: 1.15

### Turbulence Model
- Model: k-epsilon RANS
- Wall functions: Standard
- Roughness length (z0): 0.5m (urban)

### Boundary Conditions
- Inlet: Atmospheric boundary layer profile
- Outlet: Zero gradient
- Ground/Buildings: No-slip walls
- Top/Sides: Slip conditions

## Troubleshooting

### Common Issues

1. **Insufficient Memory**
   - Increase Docker memory allocation
   - Reduce mesh resolution
   - Use coarser initial mesh

2. **Convergence Problems**
   - Check mesh quality (`checkMesh`)
   - Reduce relaxation factors
   - Increase simulation steps

3. **Missing Dependencies**
   ```bash
   # In Docker container
   apt-get update
   apt-get install python3-pip
   pip3 install numpy pandas
   ```

4. **Coordinate System Issues**
   - Ensure .gpkg uses SWEREF99 TM (EPSG:3011)
   - Check domain centering in preprocessing

## Performance Optimization

### Mesh Optimization
- Use structured mesh where possible
- Apply refinement only near buildings
- Use wall functions (y+ = 30-300)

### Parallel Processing
```bash
# Decompose domain for parallel processing
decomposePar
mpirun -np 4 simpleFoam -parallel
reconstructPar
```

### Simplified Approach for Large Areas
- Run 2D simulations first
- Use periodic boundary conditions
- Implement nested domains

## Validation

### Recommended Validation Steps
1. Compare with local weather station data
2. Check mass conservation
3. Perform grid independence study
4. Validate against field measurements

### Quality Metrics
- Residuals < 10^-4
- Continuity errors < 10^-5
- y+ values in acceptable range

## References

1. OpenFOAM User Guide: https://www.openfoam.com/documentation/user-guide
2. SMHI Wind Data: https://www.smhi.se/en/services/data-and-statistics
3. Urban Wind Comfort Criteria: Lawson & Penwarden (1975)
4. CFD Best Practice Guidelines: COST Action 732

## Contact & Support

For questions about this implementation:
- Review OpenFOAM forums: https://www.cfd-online.com/Forums/openfoam/
- Check SMHI documentation for wind data
- Consult urban CFD literature for methodology

## Next Steps

After successful wind field simulation:
1. Integrate with tree canopy data (10m resolution)
2. Combine with LST (Land Surface Temperature) data
3. Apply to RAHV (Residence and Activity-based Heat Vulnerability) analysis
4. Generate urban comfort maps for planning decisions

---

*Note: This workflow is designed for research purposes. For operational use, consider professional CFD software and validation against local measurements.*
