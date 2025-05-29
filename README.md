[Paper Link](https://www.techrxiv.org/users/788558/articles/1022587-kinematic-and-dynamic-modeling-of-cable-object-interference-and-wrapping-in-complex-geometrical-shaped-cable-driven-parallel-robots)

## Installation

To set up the repository, follow these steps:

1. **Clone the repository**:
    ```bash
    git clone https://github.com/bhattner143/GeoWrapSim-CDPR.git
    ```

2. **Add the folder to MATLAB**:
    Open MATLAB, navigate to the downloaded folder, and add it to the MATLAB path.

3. **Initialize the CASPR environment**:
    Run the `initialise_CASPR.m` script to configure the Cable-Driven Parallel Robot Simulator (CASPR) environment.

### About CASPR
CASPR is a simulation framework designed for modeling and analyzing cable-driven parallel robots. It provides tools for kinematic and dynamic analysis, enabling researchers and engineers to study complex cable-driven systems effectively.
### Tutorials for CASPR

To get started with CASPR, explore the tutorial scripts provided in the `scripts/CASPR_tutorials` folder. These `.mlx` scripts offer step-by-step guidance on using CASPR for various applications, including kinematic and dynamic analysis of cable-driven parallel robots. Open these scripts in MATLAB to follow the tutorials interactively.
### Note on Tutorials

The tutorial scripts provided in the `scripts/CASPR_tutorials` folder are untested due to modifications in the repository. If you encounter any issues, you can refer to the original CASPR repository for tutorials at [CASPR GitHub Repository](https://github.com/darwinlau/CASPR).

Additionally, you can explore the CASPR tutorial series on YouTube for step-by-step guidance:  
[CASPR YouTube Playlist](https://www.youtube.com/watch?v=b_24t_j1uQo&list=PLZCfv3Lre4aVbsS8zFYlg2x-kPuB-rrRC)

## Model
## Model

### NURBS Teapot Model
The `.mat` files for the NURBS teapot model are located in the `data/nurbs_related` directory.

### BMWrapArm Model
In this work, we have utilized the BMWrapArm model. The configuration files for this model are available in the `data/model_config/models/MCDM/BMWrapArm` directory. Below is a description of the key files:

1. **BMWrapArm_bodies.xml**  
    Contains information about the CDPR bodies, including:
    - Initial joint pose
    - Center of mass and its location
    - Inertial matrix

2. **BMWrapArm_cables.xml**  
    Provides details about cable locations and attachment points.

3. **BMWrapArm_operational_spaces.xml**  
    Specifies information related to operational spaces.

4. **BMWrapArm_trajectories.xml**  
    Includes details about trajectories.

### Note on Other Robots
Other robot models are available in the `SCDM`, `MCDM`, and `HCDM` directories. However, wrapping models have not been implemented for these robots.


## Running the scripts