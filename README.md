# Lattice Boltzmann (LBM) Solver for HPC Fluid Dynamics
![C++](https://img.shields.io/badge/C%2B%2B-20-blue.svg)
![OpenMP](https://img.shields.io/badge/Parallel-OpenMP-orange.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

A parallel C++20 implementation of the Lattice Boltzmann Method optimized for High Performance Computing environments.

## Overview

This project provides a flexible framework for simulating fluid flows using the discrete Boltzmann equation. It is designed with **High Performance Computing (HPC)** principles in mind, utilizing modern C++ features (e.g. concepts, variadic templates), cache-friendly data structures, and shared-memory parallelization (OpenMP).

The solver currently supports **Lid-Driven Cavity** flow benchmarks but is structured to be extensible for other geometries and boundary conditions. 

## Key Features

* **Lattice Models:** Support for **D2Q9** (2D) and **D3Q19** (3D) lattice discrete velocity models.
* **HPC Optimizations:**
    * **Structure of Arrays (SoA)** memory layout for efficient vectorization and cache usage.
    * **OpenMP** parallelization for multi-core execution.
    * **Template classes** for maximum generality
* **Numerical Methods:**
    * **BGK** (Bhatnagar-Gross-Krook) collision operator.
    * **Bounce-Back** boundary conditions for no-slip walls.
    * **Moving Wall** boundary conditions (Lid-Driven).
* **Visualization:** **VTK (.vtk)** output support, fully compatible with [ParaView](https://www.paraview.org/).
* **Validation:** Includes tests against the **Ghia et al. (1982)** benchmark data.

## Project Structure

``` text
lbm-1-lbm/
├── include/
│   └── LBM/                      # Header-only template library (Lattice, Solver, Boundary, Descriptor, LBM)
├── scripts
|   └── job_scripts.sh            # Reference script to launch the simulation on the cluster (qsub) 
├── src/
│   └── LBM/                      # Header implementation 
|   └── main2D.cpp                # Main local simulation entry point (LBM target) in 2D
|   └── main3D.cpp                # Main local simulation entry point (LBM target) in 3D
|   └── cluster_main.cpp          # Cluster simulation entry point (LBM_cluster target) in 3D
├── test/
│   ├── lattice_test_D2Q9.cpp     # Unit tests for core data structures (Catch2)
│   └── test_ghia_Re100_2D.cpp    # Physics validation (Re=100)
├── CMakeLists.txt                # CMake build configuration
└── README.md
```

## Software Architecture

The project is designed around three main components, leveraging **C++20 Concepts** and **Static Polymorphism**.

### Class Diagram

```mermaid
classDiagram
    class isDescriptor {
        <<concept>>
        +int d (Dimension)
        +int q (Discrete Velocities)
    }

    class D2Q9 {
        +static constexpr d = 2
        +static constexpr q = 9
        +s.c. array~array[2], 9~ c
        +s.c. array~9~ w
    }

    class D3Q19 {
        +static constexpr d = 3
        +static constexpr q = 19
        +s.c. array~array[3],19~ c
        +s.c array~19~ w
    }

    class Lattice~Descriptor, float_type~ {
        +vector rho
        +array~vector~ u
        +array~vector~ f
    }

    class WallsBoundary~Descriptor, float_type~ {
        +is_moving_wall
        +wall_speed
        +is_at_bound()
        +get_speed_of_wall()
        +will_get_bounced_back()
    }

    class Solver~Descriptor, float_type~ {
        +Lattice& lattice
        +solve()
        +stream_collide()
    }

    isDescriptor <|-- D2Q9 : implements
    isDescriptor <|-- D3Q19 : implements
    Lattice ..> isDescriptor : depends on
    Lattice *-- WallsBoundary : composition
    Solver --> Lattice : manipulates
```

### Core Components

#### 1. Descriptors (`D2Q9`, `D3Q19`)
Defines the discrete velocity set ($\vec{c}_i$) and weights ($w_i$) at **compile-time**.
- Uses `constexpr` arrays to allow loop unrolling and aggressive compiler optimizations.
- Allow the solver to switch between 2D and 3D logic without runtime overhead.

#### 2. `WallsBoundary` (Geometry Handler)
Handles the geometric logic of the domain boundaries.
- **Responsibility:** Determines if a specific cell index corresponds to a physical wall and handles boundary properties (e.g., identifying the Moving Lid vs. Static Walls).

#### 3. `Lattice` (Data Container)
The owner of the simulation data and grid topology.
- **Data Layout:** Implements a **Structure of Arrays (SoA)** layout for macroscopic fields (`u`) and distributions (`f`). This layout is chosen to maximize cache locality and enable auto-vectorization (SIMD) on modern CPUs.
- **Composition:** It owns an instance of `WallsBoundary` to delegate geometric checks while maintaining a unified interface for the solver.

#### 4. `Solver` (Algorithm)
Implements the time-evolution loop of the LBM equation.
- **Fused Kernel:** Performs the *Stream* and *Collide* steps in a single pass over the grid. This maximizes arithmetic intensity and reduces memory bandwidth pressure compared to separate passes.
- **Parallelism:** The main loop is fully parallelized using **OpenMP** with a static scheduling strategy.

## Building
- C++20 or later compiler
- CMake 3.28+
- Git

### Building Instructions
``` bash
cd lbm-1-lbm
mkdir -p build
cd build
cmake .. -DBUILD_TESTS=ON
make
```

## Usage
### Testing
We use [Catch2](https://github.com/catchorg/Catch2) for unit and integration testing.
To run the tests:
``` bash
cd build/test
./run_tests
```
#### Benchmark: Ghia et al. (1982)

The project includes a specific validation test that compares the computed velocity profiles (u-velocity along vertical centerline, v-velocity along horizontal centerline) against the tabulated data from the standard reference.

### Standard usage
To maintain a clean codebase and optimize for different environments, the project provides two distinct entry points:

* **`main.cpp` (Local/Debug):** The standard entry point intended for local development, debugging, and visualization tests. Use this target for runs on personal machines or when immediate graphical feedback is required.
    
* **`cluster_main.cpp` (HPC/Production):** Specifically designed for high-performance execution on the cluster. 

This separation ensures that cluster-specific configurations do not clutter the local development logic.

## Visualization

The  generates data in **VTK format** (`.vtk`), which allows for high-quality 2D/3D visualization.

### Output Location
By default, all output files are generated **in the same directory where the executable is launched**. 
* If you run `./LBM_cluster`, the `.vtk` files will appear in the current folder.
* **Note:** For long runs on the cluster, ensure you are running the job from a directory with sufficient storage quota.

### How to Visualize (ParaView)
We recommend using [ParaView](https://www.paraview.org/) to inspect the results.

1.  **Load the Data:** Open ParaView and navigate to the build/execution directory. ParaView will automatically group the file sequence (e.g., `output_..vtk`). Select the group and click **Apply**.
2.  **Inspect Internal Flow:** If the  is 3D, the domain will appear as a solid block by default. To see the internal fluid dynamics, use standard filters:
    * **Slice:** To view 2D cross-sections of velocity or density.
    * **Stream Tracer:** To visualize flow pathlines.
##  Example results
* [2D Lid-Driven Cavity velocity vector field](https://youtu.be/lDILUFNHXN0)
* [3D Lid-Driven Cavity (Re=400)](https://youtu.be/Xj1hBF0MpeU)


## Acknowledgments
This project was developed for the **Advanced Methods for Scientific Computing** course at **Politecnico di Milano** (A.Y. 2025-2026).

* **[Prof. Luca Formaggia](https://github.com/lformaggia)** - *Course Instructor*
* **[Dr. Paolo Joseph Baioni](https://github.com/pjbaioni)** - *Teaching Assistant*
* **[Dr. Marco Scarpelli](https://github.com/ScarpMarc)** - *Teaching Assistant*

The simulation used in the presentation for this project were, in part, performed on the HPC Cluster of the Department of Mathematics of Politecnico di Milano which was funded by MUR grant Dipartimento di Eccellenza 2023-2027.

## Authors
* **[Giovanni Carpenedo](https://github.com/gcarpenedo)**
* **[Daniele Confalonieri](https://github.com/DanieleConfalonieri)**
* **[Cristiano Corona](https://github.com/CristianoCorona)**
* **[Simone Ferri](https://github.com/SimoFerri)**
* **[Federico Pizzolato](https://github.com/federico-pizz)**


