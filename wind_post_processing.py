#!/usr/bin/env python3
"""
Post-processing Script for Urban Wind Field Simulation
Extracts wind velocity at specified height and exports to GeoTIFF
"""

import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
import json
import rasterio
from rasterio.transform import from_origin
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import geopandas as gpd
from shapely.geometry import Point

class WindFieldPostProcessor:
    """
    Post-processing tools for OpenFOAM wind simulation results
    """
    
    def __init__(self, case_dir, output_dir="results"):
        self.case_dir = Path(case_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load domain information
        domain_info_file = self.case_dir.parent / "preprocessing" / "domain_info.json"
        if domain_info_file.exists():
            with open(domain_info_file, "r") as f:
                self.domain_info = json.load(f)
        else:
            self.domain_info = None
            
    def create_sampling_dict(self, height=1.5, resolution=10):
        """
        Create sampleDict for OpenFOAM to extract data at specified height
        
        Args:
            height: Height above ground to sample (meters)
            resolution: Grid resolution for sampling (meters)
        """
        print(f"Creating sampling dictionary for height={height}m, resolution={resolution}m...")
        
        if not self.domain_info:
            raise ValueError("Domain information not found. Run preprocessing first.")
        
        bounds = self.domain_info["domain_bounds"]
        center = self.domain_info["center"]
        
        # Calculate grid points
        x_min = bounds[0] - center[0]
        x_max = bounds[2] - center[0]
        y_min = bounds[1] - center[1]
        y_max = bounds[3] - center[1]
        
        nx = int((x_max - x_min) / resolution)
        ny = int((y_max - y_min) / resolution)
        
        sample_dict = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      sampleDict;
}}

type            sets;
libs            ("libsampling.so");

interpolationScheme cellPoint;

setFormat       raw;

sets
(
);

surfaces
(
    pedestrianLevel
    {{
        type            plane;
        planeType       pointAndNormal;
        pointAndNormalDict
        {{
            point   (0 0 {height});
            normal  (0 0 1);
        }}
        
        interpolate     true;
    }}
);

fields          (U p k epsilon);
"""
        
        # Alternative: Create a cutting plane dictionary for more control
        cutting_plane_dict = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      cuttingPlaneDict;
}}

type            surfaces;
libs            ("libsampling.so");

writeControl    writeTime;

surfaceFormat   vtk;

formatOptions
{{
    default
    {{
        format  ascii;
    }}
}}

interpolationScheme cellPoint;

surfaces
(
    zPlane_{height}m
    {{
        type            cuttingPlane;
        planeType       pointAndNormal;
        pointAndNormalDict
        {{
            point   (0 0 {height});
            normal  (0 0 1);
        }}
        interpolate     true;
        fields          (U p k);
    }}
);

fields          (U p k epsilon);
"""
        
        # Save sampling dictionaries
        with open(self.case_dir / "system" / "sampleDict", "w") as f:
            f.write(sample_dict)
            
        with open(self.case_dir / "system" / "cuttingPlaneDict", "w") as f:
            f.write(cutting_plane_dict)
            
        print(f"✓ Created sampling dictionaries")
        
        return nx, ny, (x_min, x_max, y_min, y_max)
    
    def extract_wind_field_paraview(self, time_step="latestTime", height=1.5):
        """
        Extract wind field using ParaView Python interface
        
        Args:
            time_step: Time step to extract (default: "latestTime")
            height: Height above ground (meters)
        """
        print(f"Extracting wind field at height={height}m using ParaView...")
        
        paraview_script = f"""
import sys
sys.path.append('/opt/ParaView/lib/python3.9/site-packages')

from paraview.simple import *
import numpy as np

# Load OpenFOAM case
case = OpenFOAMReader(FileName='{self.case_dir}/case.foam')
case.MeshRegions = ['internalMesh']
case.CellArrays = ['U', 'p', 'k', 'epsilon']

# Update pipeline
UpdatePipeline()

# Create slice at specified height
slice1 = Slice(Input=case)
slice1.SliceType = 'Plane'
slice1.SliceType.Origin = [0.0, 0.0, {height}]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# Update
UpdatePipeline()

# Save as CSV
SaveData('{self.output_dir}/wind_field_{height}m.csv', proxy=slice1)

# Save as VTK
SaveData('{self.output_dir}/wind_field_{height}m.vtk', proxy=slice1)

print("✓ Wind field extracted successfully")
"""
        
        script_path = self.output_dir / "extract_paraview.py"
        with open(script_path, "w") as f:
            f.write(paraview_script)
            
        print(f"✓ Created ParaView extraction script: {script_path}")
        
    def convert_to_geotiff(self, csv_file, output_file, resolution=10):
        """
        Convert extracted wind field data to GeoTIFF format
        
        Args:
            csv_file: Path to CSV file with wind field data
            output_file: Output GeoTIFF file path
            resolution: Grid resolution (meters)
        """
        print(f"Converting wind field to GeoTIFF...")
        
        # Read CSV data
        df = pd.read_csv(csv_file)
        
        # Extract coordinates and velocity magnitude
        if 'Points:0' in df.columns:
            x = df['Points:0'].values
            y = df['Points:1'].values
        elif 'x' in df.columns:
            x = df['x'].values
            y = df['y'].values
        else:
            raise ValueError("Cannot find coordinate columns in CSV")
        
        # Calculate velocity magnitude
        if 'U:0' in df.columns:
            u = df['U:0'].values
            v = df['U:1'].values
            w = df['U:2'].values
            velocity_mag = np.sqrt(u**2 + v**2 + w**2)
        elif 'U_mag' in df.columns:
            velocity_mag = df['U_mag'].values
        else:
            raise ValueError("Cannot find velocity columns in CSV")
        
        # Add back geographic coordinates
        if self.domain_info:
            x_geo = x + self.domain_info["center"][0]
            y_geo = y + self.domain_info["center"][1]
        else:
            x_geo = x
            y_geo = y
        
        # Create regular grid
        x_min, x_max = x_geo.min(), x_geo.max()
        y_min, y_max = y_geo.min(), y_geo.max()
        
        nx = int((x_max - x_min) / resolution)
        ny = int((y_max - y_min) / resolution)
        
        xi = np.linspace(x_min, x_max, nx)
        yi = np.linspace(y_min, y_max, ny)
        xi_grid, yi_grid = np.meshgrid(xi, yi)
        
        # Interpolate to regular grid
        zi = griddata(
            (x_geo, y_geo), 
            velocity_mag, 
            (xi_grid, yi_grid), 
            method='linear',
            fill_value=0
        )
        
        # Create GeoTIFF
        transform = from_origin(x_min, y_max, resolution, resolution)
        
        with rasterio.open(
            output_file,
            'w',
            driver='GTiff',
            height=ny,
            width=nx,
            count=1,
            dtype=zi.dtype,
            crs='EPSG:3011',  # SWEREF99 TM for Stockholm
            transform=transform
        ) as dst:
            dst.write(zi, 1)
        
        print(f"✓ Saved GeoTIFF: {output_file}")
        
        # Also save velocity components if needed
        if 'U:0' in df.columns:
            # Save U component
            ui = griddata((x_geo, y_geo), u, (xi_grid, yi_grid), method='linear', fill_value=0)
            u_file = str(output_file).replace('.tif', '_u.tif')
            with rasterio.open(
                u_file, 'w', driver='GTiff', height=ny, width=nx,
                count=1, dtype=ui.dtype, crs='EPSG:3011', transform=transform
            ) as dst:
                dst.write(ui, 1)
            
            # Save V component
            vi = griddata((x_geo, y_geo), v, (xi_grid, yi_grid), method='linear', fill_value=0)
            v_file = str(output_file).replace('.tif', '_v.tif')
            with rasterio.open(
                v_file, 'w', driver='GTiff', height=ny, width=nx,
                count=1, dtype=vi.dtype, crs='EPSG:3011', transform=transform
            ) as dst:
                dst.write(vi, 1)
                
            print(f"✓ Saved velocity component GeoTIFFs")
        
        return zi, transform
    
    def combine_wind_directions(self, wind_results_dict, output_file):
        """
        Combine wind field results from multiple directions based on probability weighting
        
        Args:
            wind_results_dict: Dictionary with wind direction results and probabilities
            output_file: Output combined GeoTIFF file
        """
        print("Combining wind field results from multiple directions...")
        
        combined_field = None
        total_weight = 0
        
        for direction, data in wind_results_dict.items():
            tif_file = self.output_dir / f"wind_{direction}_{data['angle']}deg.tif"
            
            if tif_file.exists():
                with rasterio.open(tif_file) as src:
                    wind_field = src.read(1)
                    weight = data['probability']
                    
                    if combined_field is None:
                        combined_field = np.zeros_like(wind_field)
                        transform = src.transform
                        profile = src.profile
                    
                    combined_field += wind_field * weight
                    total_weight += weight
        
        # Normalize by total weight
        if total_weight > 0:
            combined_field /= total_weight
        
        # Save combined field
        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(combined_field, 1)
        
        print(f"✓ Saved combined wind field: {output_file}")
        
        return combined_field
    
    def create_visualization(self, geotiff_file, output_image):
        """
        Create visualization of wind field
        
        Args:
            geotiff_file: Path to GeoTIFF file
            output_image: Output image file path
        """
        print("Creating wind field visualization...")
        
        # Read GeoTIFF
        with rasterio.open(geotiff_file) as src:
            wind_field = src.read(1)
            transform = src.transform
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Create custom colormap
        colors = ['#0000ff', '#00ffff', '#00ff00', '#ffff00', '#ff0000']
        n_bins = 100
        cmap = LinearSegmentedColormap.from_list('wind_speed', colors, N=n_bins)
        
        # Plot wind field
        im = ax.imshow(
            wind_field,
            cmap=cmap,
            extent=[
                transform[2], transform[2] + transform[0] * wind_field.shape[1],
                transform[5] + transform[4] * wind_field.shape[0], transform[5]
            ],
            origin='upper',
            interpolation='bilinear'
        )
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, label='Wind Speed (m/s)')
        
        # Labels and title
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        ax.set_title('Urban Wind Field at Pedestrian Level (1.5m)')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        stats_text = f"Min: {wind_field.min():.2f} m/s\n"
        stats_text += f"Max: {wind_field.max():.2f} m/s\n"
        stats_text += f"Mean: {wind_field.mean():.2f} m/s\n"
        stats_text += f"Std: {wind_field.std():.2f} m/s"
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                verticalalignment='top', bbox=dict(boxstyle='round', 
                facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Saved visualization: {output_image}")
    
    def export_comfort_zones(self, wind_field_file, output_shapefile):
        """
        Export wind comfort zones as shapefile based on Lawson criteria
        
        Args:
            wind_field_file: Path to wind field GeoTIFF
            output_shapefile: Output shapefile path
        """
        print("Calculating wind comfort zones...")
        
        # Lawson comfort criteria (m/s)
        comfort_criteria = {
            "sitting": 4.0,
            "standing": 6.0,
            "walking": 8.0,
            "uncomfortable": 10.0,
            "dangerous": 15.0
        }
        
        # Read wind field
        with rasterio.open(wind_field_file) as src:
            wind_field = src.read(1)
            transform = src.transform
        
        # Create comfort zone classification
        comfort_zones = np.zeros_like(wind_field, dtype=int)
        comfort_zones[wind_field <= comfort_criteria["sitting"]] = 1  # Sitting
        comfort_zones[(wind_field > comfort_criteria["sitting"]) & 
                     (wind_field <= comfort_criteria["standing"])] = 2  # Standing
        comfort_zones[(wind_field > comfort_criteria["standing"]) & 
                     (wind_field <= comfort_criteria["walking"])] = 3  # Walking
        comfort_zones[(wind_field > comfort_criteria["walking"]) & 
                     (wind_field <= comfort_criteria["uncomfortable"])] = 4  # Uncomfortable
        comfort_zones[wind_field > comfort_criteria["uncomfortable"]] = 5  # Dangerous
        
        # Save as new GeoTIFF
        comfort_tiff = str(output_shapefile).replace('.shp', '_zones.tif')
        profile = src.profile
        profile.update(dtype=rasterio.uint8, count=1)
        
        with rasterio.open(comfort_tiff, 'w', **profile) as dst:
            dst.write(comfort_zones.astype(np.uint8), 1)
        
        print(f"✓ Saved comfort zones GeoTIFF: {comfort_tiff}")
        
        # Calculate statistics
        unique, counts = np.unique(comfort_zones, return_counts=True)
        total_pixels = comfort_zones.size
        
        print("\nWind Comfort Statistics:")
        print("-" * 40)
        comfort_names = ["", "Sitting", "Standing", "Walking", "Uncomfortable", "Dangerous"]
        for zone, count in zip(unique, counts):
            if zone > 0:
                percentage = (count / total_pixels) * 100
                print(f"{comfort_names[zone]:15} {percentage:6.2f}%")
        
        return comfort_zones
    
    def run_post_processing(self, height=1.5, resolution=10):
        """
        Run complete post-processing workflow
        
        Args:
            height: Height above ground for extraction (meters)
            resolution: Grid resolution (meters)
        """
        print("\n" + "="*60)
        print("POST-PROCESSING WORKFLOW")
        print("="*60 + "\n")
        
        # Create sampling dictionary
        nx, ny, bounds = self.create_sampling_dict(height, resolution)
        
        # Extract using ParaView
        self.extract_wind_field_paraview(height=height)
        
        # Note: The actual extraction would be done using OpenFOAM tools
        print("\nTo extract data using OpenFOAM:")
        print(f"1. cd {self.case_dir}")
        print("2. postProcess -func sampleDict")
        print("3. postProcess -func cuttingPlaneDict")
        
        # Create example workflow for GeoTIFF conversion
        csv_example = self.output_dir / f"wind_field_{height}m.csv"
        tiff_output = self.output_dir / f"wind_field_{height}m.tif"
        
        print(f"\nAfter extraction, convert to GeoTIFF:")
        print(f"python -c \"from wind_post_processing import WindFieldPostProcessor")
        print(f"processor = WindFieldPostProcessor('{self.case_dir}', '{self.output_dir}')")
        print(f"processor.convert_to_geotiff('{csv_example}', '{tiff_output}', {resolution})\"")
        
        # Create visualization example
        image_output = self.output_dir / f"wind_field_{height}m.png"
        print(f"\nCreate visualization:")
        print(f"processor.create_visualization('{tiff_output}', '{image_output}')")
        
        # Export comfort zones
        shapefile_output = self.output_dir / f"comfort_zones_{height}m.shp"
        print(f"\nExport comfort zones:")
        print(f"processor.export_comfort_zones('{tiff_output}', '{shapefile_output}')")
        
        print("\n" + "="*60)
        print("POST-PROCESSING SETUP COMPLETE")
        print("="*60)

# Example usage
if __name__ == "__main__":
    # Initialize post-processor
    processor = WindFieldPostProcessor(
        case_dir="stockholm_wind_sim/openfoam_case",
        output_dir="stockholm_wind_sim/results"
    )
    
    # Run post-processing workflow
    processor.run_post_processing(height=1.5, resolution=10)
    
    # Example of combining multiple wind directions
    wind_results = {
        "W": {"angle": 270, "probability": 0.224},
        "S": {"angle": 180, "probability": 0.174},
        "SW": {"angle": 225, "probability": 0.134}
    }
    
    # processor.combine_wind_directions(
    #     wind_results, 
    #     processor.output_dir / "combined_wind_field.tif"
    # )
