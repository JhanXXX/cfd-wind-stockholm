#!/usr/bin/env python3
"""
Urban Wind Field Simulation Workflow for Stockholm
Author: Urban Climate Analysis Project
Description: Complete workflow for CFD wind simulation using OpenFOAM
"""

import os
import numpy as np
import geopandas as gpd
from pathlib import Path
import subprocess
import json
from datetime import datetime

class UrbanWindSimulation:
    """
    Main class for managing urban wind field CFD simulation
    """
    
    def __init__(self, project_dir="stockholm_wind_sim"):
        self.project_dir = Path(project_dir)
        self.case_dir = self.project_dir / "openfoam_case"
        
        # Stockholm wind statistics (summer months)
        self.wind_data = {
            "summer": {
                "directions": {
                    "W": {"probability": 0.224, "angle": 270, "speed": 4.0},
                    "S": {"probability": 0.174, "angle": 180, "speed": 3.5},
                    "SW": {"probability": 0.134, "angle": 225, "speed": 3.8},
                    "N": {"probability": 0.113, "angle": 0, "speed": 3.2},
                    "SE": {"probability": 0.104, "angle": 135, "speed": 3.3},
                    "NW": {"probability": 0.099, "angle": 315, "speed": 3.5},
                    "E": {"probability": 0.091, "angle": 90, "speed": 3.0},
                    "NE": {"probability": 0.062, "angle": 45, "speed": 2.8}
                }
            }
        }
        
        # Urban roughness parameters
        self.roughness_params = {
            "z0": 0.5,  # Roughness length for urban area (m)
            "d": 5.0,    # Displacement height (m)
            "zref": 10.0, # Reference height (m)
            "kappa": 0.41 # von Karman constant
        }
        
        # Simulation parameters
        self.sim_params = {
            "domain_height": 300,  # meters
            "upstream": 500,       # meters
            "downstream": 1500,    # meters
            "lateral": 500,        # meters
            "first_cell_height": 0.5,  # meters
            "expansion_ratio": 1.15,
            "target_resolution": 10  # meters
        }
        
    def setup_directories(self):
        """Create necessary directory structure"""
        dirs = [
            self.case_dir,
            self.case_dir / "0",
            self.case_dir / "constant",
            self.case_dir / "constant" / "triSurface",
            self.case_dir / "system",
            self.project_dir / "preprocessing",
            self.project_dir / "results",
            self.project_dir / "postprocessing"
        ]
        
        for d in dirs:
            d.mkdir(parents=True, exist_ok=True)
            
        print(f"✓ Directory structure created at {self.project_dir}")
        
    def process_building_data(self, gpkg_file):
        """
        Process building data from GeoPackage file
        
        Args:
            gpkg_file: Path to .gpkg file containing building polygons and heights
        """
        print("Processing building data...")
        
        # Read building data
        buildings = gpd.read_file(gpkg_file)
        
        # Assuming columns: 'geometry', 'height' 
        # Adjust column names as needed
        if 'height' not in buildings.columns:
            print("Warning: 'height' column not found. Using default height of 15m")
            buildings['height'] = 15.0
        
        # Get bounds for domain setup
        bounds = buildings.total_bounds  # [minx, miny, maxx, maxy]
        
        # Calculate domain extents
        domain_info = {
            "building_bounds": bounds.tolist(),
            "domain_bounds": [
                bounds[0] - self.sim_params["upstream"],
                bounds[1] - self.sim_params["lateral"],
                bounds[2] + self.sim_params["downstream"],
                bounds[3] + self.sim_params["lateral"]
            ],
            "center": [
                (bounds[0] + bounds[2]) / 2,
                (bounds[1] + bounds[3]) / 2
            ]
        }
        
        # Save domain info
        with open(self.project_dir / "preprocessing" / "domain_info.json", "w") as f:
            json.dump(domain_info, f, indent=2)
        
        # Export buildings to STL format
        self.export_buildings_to_stl(buildings, domain_info["center"])
        
        return domain_info
    
    def export_buildings_to_stl(self, buildings, center):
        """
        Convert building polygons to STL format for OpenFOAM
        """
        from shapely.geometry import Polygon
        import trimesh
        
        print("Converting buildings to STL...")
        
        # Translate coordinates to center origin
        buildings = buildings.copy()
        buildings['geometry'] = buildings['geometry'].translate(
            xoff=-center[0], yoff=-center[1]
        )
        
        meshes = []
        
        for idx, building in buildings.iterrows():
            if building.geometry.geom_type == 'Polygon':
                # Create 3D mesh from polygon
                height = building['height']
                polygon = building.geometry
                
                # Get exterior coordinates
                coords = list(polygon.exterior.coords)[:-1]  # Remove duplicate last point
                
                # Create vertices for bottom and top
                vertices_bottom = [[x, y, 0] for x, y in coords]
                vertices_top = [[x, y, height] for x, y in coords]
                
                # Create faces
                n = len(coords)
                faces = []
                
                # Bottom face
                faces.extend(self.triangulate_polygon(vertices_bottom, reverse=True))
                
                # Top face  
                faces.extend(self.triangulate_polygon(vertices_top, reverse=False))
                
                # Side faces
                for i in range(n):
                    j = (i + 1) % n
                    # Create two triangles for each rectangular side
                    faces.append([i, j, i + n])
                    faces.append([j, j + n, i + n])
                
                # Create mesh
                all_vertices = vertices_bottom + vertices_top
                mesh = trimesh.Trimesh(vertices=all_vertices, faces=faces)
                meshes.append(mesh)
        
        # Combine all building meshes
        combined_mesh = trimesh.util.concatenate(meshes)
        
        # Export to STL
        stl_path = self.case_dir / "constant" / "triSurface" / "buildings.stl"
        combined_mesh.export(str(stl_path))
        
        print(f"✓ Exported {len(buildings)} buildings to STL")
        
    def triangulate_polygon(self, vertices, reverse=False):
        """Simple ear clipping triangulation for convex polygons"""
        n = len(vertices)
        if n < 3:
            return []
        
        faces = []
        for i in range(1, n - 1):
            if reverse:
                faces.append([0, i + 1, i])
            else:
                faces.append([0, i, i + 1])
        
        return faces
    
    def create_blockmesh_dict(self, domain_info):
        """
        Create blockMeshDict for OpenFOAM
        """
        print("Creating blockMeshDict...")
        
        bounds = domain_info["domain_bounds"]
        height = self.sim_params["domain_height"]
        
        # Calculate cell numbers
        dx = self.sim_params["target_resolution"]
        nx = int((bounds[2] - bounds[0]) / dx)
        ny = int((bounds[3] - bounds[1]) / dx)
        
        # Vertical cells with expansion
        nz = 50  # Number of vertical cells
        
        blockmesh_dict = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    object      blockMeshDict;
}}

scale   1;

vertices
(
    ({bounds[0] - domain_info['center'][0]} {bounds[1] - domain_info['center'][1]} 0)
    ({bounds[2] - domain_info['center'][0]} {bounds[1] - domain_info['center'][1]} 0)
    ({bounds[2] - domain_info['center'][0]} {bounds[3] - domain_info['center'][1]} 0)
    ({bounds[0] - domain_info['center'][0]} {bounds[3] - domain_info['center'][1]} 0)
    ({bounds[0] - domain_info['center'][0]} {bounds[1] - domain_info['center'][1]} {height})
    ({bounds[2] - domain_info['center'][0]} {bounds[1] - domain_info['center'][1]} {height})
    ({bounds[2] - domain_info['center'][0]} {bounds[3] - domain_info['center'][1]} {height})
    ({bounds[0] - domain_info['center'][0]} {bounds[3] - domain_info['center'][1]} {height})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) 
    simpleGrading (1 1 {self.sim_params['expansion_ratio']})
);

edges
(
);

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }}
    
    outlet
    {{
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }}
    
    ground
    {{
        type wall;
        faces
        (
            (0 1 2 3)
        );
    }}
    
    top
    {{
        type patch;
        faces
        (
            (4 5 6 7)
        );
    }}
    
    sides
    {{
        type patch;
        faces
        (
            (0 1 5 4)
            (2 3 7 6)
        );
    }}
);

mergePatchPairs
(
);
"""
        
        # Save blockMeshDict
        with open(self.case_dir / "system" / "blockMeshDict", "w") as f:
            f.write(blockmesh_dict)
            
        print(f"✓ Created blockMeshDict with {nx}x{ny}x{nz} cells")
        
    def create_snappyhexmesh_dict(self):
        """
        Create snappyHexMeshDict for OpenFOAM
        """
        print("Creating snappyHexMeshDict...")
        
        snappy_dict = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}

castellatedMesh true;
snap            true;
addLayers       true;

geometry
{
    buildings.stl
    {
        type triSurfaceMesh;
        name buildings;
    }
}

castellatedMeshControls
{
    maxLocalCells 10000000;
    maxGlobalCells 20000000;
    minRefinementCells 10;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;

    features
    (
        {
            file "buildings.eMesh";
            level 2;
        }
    );

    refinementSurfaces
    {
        buildings
        {
            level (2 3);
            patchInfo
            {
                type wall;
            }
        }
    }

    resolveFeatureAngle 30;

    refinementRegions
    {
        buildings
        {
            mode distance;
            levels ((50 2) (100 1));
        }
    }

    locationInMesh (0 0 50);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 100;
    nRelaxIter 5;
    nFeatureSnapIter 10;
    implicitFeatureSnap false;
    explicitFeatureSnap true;
    multiRegionFeatureSnap false;
}

addLayersControls
{
    relativeSizes true;

    layers
    {
        buildings
        {
            nSurfaceLayers 3;
        }
    }

    expansionRatio 1.2;
    finalLayerThickness 0.5;
    minThickness 0.25;
    nGrow 0;

    featureAngle 60;
    slipFeatureAngle 30;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
    nRelaxedIter 20;
}

meshQualityControls
{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality 1e-15;
    minArea -1;
    minTwist 0.05;
    minDeterminant 0.001;
    minFaceWeight 0.05;
    minVolRatio 0.01;
    minTriangleTwist -1;

    nSmoothScale 4;
    errorReduction 0.75;

    relaxed
    {
        maxNonOrtho 75;
    }
}

writeFlags
(
    scalarLevels
    layerSets
    layerFields
);

mergeTolerance 1e-6;
"""
        
        with open(self.case_dir / "system" / "snappyHexMeshDict", "w") as f:
            f.write(snappy_dict)
            
        print("✓ Created snappyHexMeshDict")

    def create_boundary_conditions(self, wind_direction, wind_speed):
        """
        Create boundary condition files for OpenFOAM
        
        Args:
            wind_direction: Wind direction in degrees (0=North, 90=East, etc.)
            wind_speed: Reference wind speed at 10m height (m/s)
        """
        print(f"Creating boundary conditions for {wind_direction}° at {wind_speed} m/s...")
        
        # Convert wind direction to velocity components
        # OpenFOAM uses "wind coming from" convention
        angle_rad = np.radians(wind_direction + 180)  # Add 180 to get "wind going to"
        u = wind_speed * np.cos(angle_rad)
        v = wind_speed * np.sin(angle_rad)
        
        # Create U (velocity) file
        u_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({u} {v} 0);

boundaryField
{{
    inlet
    {{
        type            atmBoundaryLayerInletVelocity;
        Uref            {wind_speed};
        Zref            {self.roughness_params['zref']};
        zDir            (0 0 1);
        flowDir         ({np.cos(angle_rad)} {np.sin(angle_rad)} 0);
        z0              uniform {self.roughness_params['z0']};
        zGround         uniform 0.0;
        value           uniform ({u} {v} 0);
    }}

    outlet
    {{
        type            zeroGradient;
    }}

    ground
    {{
        type            noSlip;
    }}

    buildings
    {{
        type            noSlip;
    }}

    top
    {{
        type            slip;
    }}

    sides
    {{
        type            slip;
    }}
}}
"""
        
        # Create p (pressure) file
        p_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }

    outlet
    {
        type            fixedValue;
        value           uniform 0;
    }

    ground
    {
        type            zeroGradient;
    }

    buildings
    {
        type            zeroGradient;
    }

    top
    {
        type            zeroGradient;
    }

    sides
    {
        type            zeroGradient;
    }
}
"""
        
        # Create k (turbulent kinetic energy) file
        k_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       volScalarField;
    object      k;
}}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0.1;

boundaryField
{{
    inlet
    {{
        type            atmBoundaryLayerInletK;
        z0              uniform {self.roughness_params['z0']};
        Uref            {wind_speed};
        Zref            {self.roughness_params['zref']};
        zDir            (0 0 1);
        flowDir         ({np.cos(angle_rad)} {np.sin(angle_rad)} 0);
        zGround         uniform 0.0;
        value           uniform 0.1;
    }}

    outlet
    {{
        type            zeroGradient;
    }}

    ground
    {{
        type            kqRWallFunction;
        value           uniform 0.1;
    }}

    buildings
    {{
        type            kqRWallFunction;
        value           uniform 0.1;
    }}

    top
    {{
        type            slip;
    }}

    sides
    {{
        type            slip;
    }}
}}
"""
        
        # Create epsilon (turbulent dissipation) file
        epsilon_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       volScalarField;
    object      epsilon;
}}

dimensions      [0 2 -3 0 0 0 0];

internalField   uniform 0.01;

boundaryField
{{
    inlet
    {{
        type            atmBoundaryLayerInletEpsilon;
        z0              uniform {self.roughness_params['z0']};
        Uref            {wind_speed};
        Zref            {self.roughness_params['zref']};
        zDir            (0 0 1);
        flowDir         ({np.cos(angle_rad)} {np.sin(angle_rad)} 0);
        zGround         uniform 0.0;
        value           uniform 0.01;
    }}

    outlet
    {{
        type            zeroGradient;
    }}

    ground
    {{
        type            epsilonWallFunction;
        value           uniform 0.01;
    }}

    buildings
    {{
        type            epsilonWallFunction;
        value           uniform 0.01;
    }}

    top
    {{
        type            slip;
    }}

    sides
    {{
        type            slip;
    }}
}}
"""
        
        # Create nut (turbulent viscosity) file
        nut_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            calculated;
        value           uniform 0;
    }

    outlet
    {
        type            calculated;
        value           uniform 0;
    }

    ground
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    buildings
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    top
    {
        type            calculated;
        value           uniform 0;
    }

    sides
    {
        type            calculated;
        value           uniform 0;
    }
}
"""
        
        # Save all boundary condition files
        bc_files = {
            "U": u_file,
            "p": p_file,
            "k": k_file,
            "epsilon": epsilon_file,
            "nut": nut_file
        }
        
        for filename, content in bc_files.items():
            with open(self.case_dir / "0" / filename, "w") as f:
                f.write(content)
        
        print(f"✓ Created boundary conditions for wind from {wind_direction}°")

    def create_control_dicts(self):
        """
        Create control dictionaries for OpenFOAM
        """
        print("Creating control dictionaries...")
        
        # controlDict
        control_dict = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application     simpleFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         2000;

deltaT          1;

writeControl    timeStep;

writeInterval   100;

purgeWrite      2;

writeFormat     binary;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

functions
{
    residuals
    {
        type            residuals;
        libs            ("libutilityFunctionObjects.so");
        writeControl    timeStep;
        writeInterval   1;
        fields          (U p k epsilon);
    }
}
"""
        
        # fvSchemes
        fv_schemes = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss limitedLinear 1;
    div(phi,epsilon) bounded Gauss limitedLinear 1;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}
"""
        
        # fvSolution
        fv_solution = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.1;
        smoother        GaussSeidel;
    }

    "(U|k|epsilon)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;

    residualControl
    {
        p               1e-4;
        U               1e-4;
        "(k|epsilon)"   1e-4;
    }
}

relaxationFactors
{
    equations
    {
        U               0.7;
        k               0.7;
        epsilon         0.7;
    }
}
"""
        
        # transportProperties
        transport_properties = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}

transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] 1.5e-05;
"""
        
        # turbulenceProperties
        turbulence_properties = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}

simulationType  RAS;

RAS
{
    RASModel        kEpsilon;

    turbulence      on;

    printCoeffs     on;
}
"""
        
        # Save all control dictionaries
        with open(self.case_dir / "system" / "controlDict", "w") as f:
            f.write(control_dict)
        
        with open(self.case_dir / "system" / "fvSchemes", "w") as f:
            f.write(fv_schemes)
        
        with open(self.case_dir / "system" / "fvSolution", "w") as f:
            f.write(fv_solution)
        
        with open(self.case_dir / "constant" / "transportProperties", "w") as f:
            f.write(transport_properties)
        
        with open(self.case_dir / "constant" / "turbulenceProperties", "w") as f:
            f.write(turbulence_properties)
        
        print("✓ Created all control dictionaries")

    def create_run_script(self):
        """
        Create bash script to run OpenFOAM simulation
        """
        run_script = """#!/bin/bash
#
# OpenFOAM Urban Wind Simulation Run Script
# 

echo "Starting Urban Wind Field Simulation..."
echo "======================================"

# Source OpenFOAM bashrc
source /opt/openfoam11/etc/bashrc

# Change to case directory
cd """ + str(self.case_dir) + """

# Clean previous results (optional)
if [ "$1" == "clean" ]; then
    echo "Cleaning previous results..."
    foamListTimes -rm
    rm -rf constant/polyMesh
    rm -rf processor*
    rm -f log.*
fi

# Create base mesh
echo "Step 1: Creating base mesh with blockMesh..."
blockMesh > log.blockMesh 2>&1
if [ $? -ne 0 ]; then
    echo "Error in blockMesh. Check log.blockMesh"
    exit 1
fi

# Extract surface features
echo "Step 2: Extracting surface features..."
surfaceFeatureExtract > log.surfaceFeatureExtract 2>&1

# Refine mesh around buildings
echo "Step 3: Refining mesh with snappyHexMesh..."
snappyHexMesh -overwrite > log.snappyHexMesh 2>&1
if [ $? -ne 0 ]; then
    echo "Error in snappyHexMesh. Check log.snappyHexMesh"
    exit 1
fi

# Check mesh quality
echo "Step 4: Checking mesh quality..."
checkMesh > log.checkMesh 2>&1

# Run simulation
echo "Step 5: Running simpleFoam solver..."
simpleFoam > log.simpleFoam 2>&1
if [ $? -ne 0 ]; then
    echo "Error in simpleFoam. Check log.simpleFoam"
    exit 1
fi

echo "======================================"
echo "Simulation completed successfully!"
echo "Results are in the case directory"
"""
        
        script_path = self.project_dir / "run_simulation.sh"
        with open(script_path, "w") as f:
            f.write(run_script)
        
        # Make script executable
        os.chmod(script_path, 0o755)
        
        print(f"✓ Created run script: {script_path}")

    def run_simulation(self, gpkg_file):
        """
        Main method to run the complete simulation workflow
        
        Args:
            gpkg_file: Path to GeoPackage file with building data
        """
        print("\n" + "="*60)
        print("URBAN WIND FIELD SIMULATION WORKFLOW")
        print("="*60 + "\n")
        
        # Setup
        self.setup_directories()
        
        # Process building data
        domain_info = self.process_building_data(gpkg_file)
        
        # Create mesh files
        self.create_blockmesh_dict(domain_info)
        self.create_snappyhexmesh_dict()
        
        # Create control dictionaries
        self.create_control_dicts()
        
        # Run simulations for main wind directions
        results = {}
        for direction, data in self.wind_data["summer"]["directions"].items():
            if data["probability"] > 0.1:  # Only simulate significant directions
                print(f"\n--- Simulating wind from {direction} ({data['angle']}°) ---")
                self.create_boundary_conditions(data["angle"], data["speed"])
                
                # Note: Actual simulation would be run here
                # subprocess.run(["bash", str(self.project_dir / "run_simulation.sh")])
                
                results[direction] = {
                    "angle": data["angle"],
                    "speed": data["speed"],
                    "probability": data["probability"]
                }
        
        # Create run script for manual execution
        self.create_run_script()
        
        print("\n" + "="*60)
        print("WORKFLOW SETUP COMPLETE")
        print("="*60)
        print(f"\nProject directory: {self.project_dir}")
        print(f"OpenFOAM case: {self.case_dir}")
        print("\nTo run simulation in Docker:")
        print("1. docker run -it -v $(pwd):/home/openfoam openfoam/openfoam11-paraview510")
        print("2. cd /home/openfoam/" + str(self.project_dir))
        print("3. bash run_simulation.sh")
        
        return results

# Example usage
if __name__ == "__main__":
    # Create simulation instance
    sim = UrbanWindSimulation("stockholm_wind_sim")
    
    # Run workflow (replace with actual .gpkg file path)
    # results = sim.run_simulation("path/to/stockholm_buildings.gpkg")
    
    print("\nNote: Replace 'path/to/stockholm_buildings.gpkg' with actual file path")
