# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

LatticeLab.jl is a Julia package for modeling lattice systems in condensed matter physics. It provides tools for creating crystal lattices, calculating band structures, dynamical matrices, and analyzing physical properties.

## Common Commands

### Package Development
- `using LatticeLab` - Load the package in Julia REPL
- `] activate .` - Activate the package environment
- `] test` - Run tests (if available)
- `] build` - Build the package

### Running Examples
Examples are located in the `/examples/` directory:
```julia
julia examples/example\ 1\ -\ create\ lattice.jl
julia examples/example\ 2\ -\ create\ lattice.jl
```

## Architecture

### Core Data Structures

1. **UnitCell** - Defines the fundamental repeating unit of a crystal lattice
   - Bravais lattice vectors
   - Sublattice positions
   - Atomic masses and orbitals

2. **Lattice** - The main data structure containing:
   - `R0`: Equilibrium positions of lattice sites
   - `BBOX`: Bounding box defining the piece of material
   - `LN`: Link information (connections between sites)
   - `f`: Connectivity matrix
   - `UC`: Unit cell data

3. **LinkInfo** - Describes connections/hopping between lattice sites

### Key Modules

- **lattice/**: Core lattice construction and manipulation
  - `build_lattice.jl`: Main lattice building functions
  - `Lattice.jl`: Main Lattice struct definition
  - `UnitCell.jl`: Unit cell definitions

- **hopping/**: Electronic Hamiltonian models
  - `HoppingHamiltonianQ.jl`: k-space tight-binding Hamiltonians
  - `HoppingHamiltonianR.jl`: Real-space Hamiltonians
  - `BdGHamiltonianQ.jl`: Bogoliubov-de Gennes Hamiltonians

- **dynmat/**: Dynamical matrices for phonons
  - `DynamicalMatrixQ.jl`: k-space dynamical matrices
  - `DynamicalMatrixR.jl`: Real-space dynamical matrices
  - `Phonopy.jl`: Interface to Phonopy force constants

- **band/**: Band structure calculations
  - `BandStructure.jl`: Band structure data structures
  - `bands_with_markers.jl`: Band plotting with markers
  - `compute_band_markers.jl`: Band character analysis

### Typical Workflow

1. **Create a UnitCell**: Define the crystal structure
   ```julia
   UC = UnitCell(dim, nsubl, bravais_basis, sublattice_coords, masses, orbits)
   ```

2. **Define LinkInfo**: Specify hopping/interaction patterns
   ```julia
   LN = link_info_by_distance_direction(distance_direction_dict, UC)
   ```

3. **Build Lattice**: Create the finite lattice structure
   ```julia
   lattice = build_lattice(LN, bounding_box)
   ```

4. **Create Hamiltonian**: Build electronic or phononic Hamiltonian
   ```julia
   HQ = kspace_hopping_hamiltonian(lattice, hopping_params)
   ```

5. **Calculate Band Structure**: Compute and analyze bands
   ```julia
   bands = band_structure(HQ, kpath)
   ```

## Dependencies

- `LinearAlgebra`: Matrix operations
- `SparseArrays`: Sparse matrix support
- `Spglib`: Space group operations
- `DataStructures`: Additional data structures
- `PlotSVG`: SVG plotting capabilities

## File Organization

- `/src/LatticeLab.jl`: Main module file with all exports
- `/src/lattice/`: Core lattice functionality
- `/src/hopping/`: Electronic structure
- `/src/dynmat/`: Phonon calculations
- `/src/band/`: Band structure analysis
- `/examples/`: Usage examples
- `/test/`: Test files