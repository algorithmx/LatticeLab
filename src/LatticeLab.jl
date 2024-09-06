__precompile__()


module LatticeLab

global const _DEBUG_MODE_ = false

using LinearAlgebra
using SparseArrays
#using JLD2
#using FileIO

using PlotSVG

using Spglib
using DataStructures

import SparseArrays.findnz

import Base: copy, convert, isequal, +, *, ==

import Base: print

export copy, convert, isequal, +, *, ==, findnz

export ⇒
export list_of_pairs_to_dict_merge_v
include("util/uniquen.jl")
include("util/tuplejoin.jl")
include("util/oneliners.jl")
include("util/dict_tools.jl")


include("lattice/error.jl")


export LatticeInfo, Var, Coordinates, Orbits
export Index, Indices, Mode, Modes, Masses
export check_compat
include("lattice/BasicTypes.jl")

export rotate
include("util/rotate.jl")

export BoundingBox, inbbox, bounding_box_Nmin_Nmax
include("lattice/BoundingBox.jl")

export position_index, check_valid_EqV
include("lattice/EqV.jl")

#+ TODO add doc string
export UnitCell, is_UnitCell
export generate_blocks_for_DQ, generate_blocks_for_HQ
export orbit_numbers, orbit_number, total_num_orbits
export reciprocal_basis
export num_sublattice, mass_sublattice
include("lattice/UnitCell.jl")

#+ TODO add doc string
export LinkInfo, Spring, Direction
include("lattice/LinkInfo.jl")

export link_info_by_distance_direction
export link_info_sublattice_pairs_by_nth
export link_info_element_pairs_by_nth
export link_info_sublattice_pairs_by_distance
export link_info_all_different
export pick_nearest_by_ij, pick_dist_by_ij
export all_links_between_sublattices
export all_links_dict_to_sp_no_directions
include("lattice/LinkInfo_constructors.jl")


#: NO CHANGE
export MillerIndex, to_projector
include("lattice/MillerIndex.jl")

export Lattice, is_Lattice
export lattice_dimensions, dimensions
export inner_site_id, inner_site_sublattices, inner_site_coords, inner_site_atom_labels
export num_sites, num_sublattice, num_inner_sites, inner_site_sublattices
export outer_border_site_id, inner_equivalent_of_outer_border_site_id
export orbit_numbers, orbit_number, total_num_orbits
export mass_site, is_inner_site, is_outer_site_at_boundary, is_outer_site_irrelavent
export equiv_site, num_unitcells
export find_index, find_index_nocheck
export block_ranges
export verify_monolayer
export reciprocal_basis, relative_reciprocal_vector
export unique_link_types
include("lattice/Lattice.jl")


#: NEW
export MinimalLattice, is_valid
export make_featured_vertex_list, make_featured_edge_list
include("lattice/MinimalLattice.jl")

#: NO CHANGE
export check_pbc
include("lattice/PBC.jl")

#: NO CHANGE
export translate_dN_units, TX, TY, TZ
include("lattice/translate.jl")

#export index_R0
#include("util/index_R0.jl")

#+ TODO
export compute_Bravais_cutoff
export index_t, iter_t
export reorganize_SPNB, build_lattice, fold_f_matrix
include("lattice/build_lattice.jl")

#+ TODO
export lattice_as_UC_FC, enlarge
include("lattice/enlarge.jl")

#: DONE
#export save_Lattice,            load_Lattice
#export save_MinimalLattice,     load_MinimalLattice 
#export save_MinimalLattice_VL_EL, save_Lattice_as_Minimal_VL_EL
#include("lattice/save_load.jl")

#: NO CHANGE
export DR_mode_to_ΔR, DQ_mode_to_ΔR
export show_lattice_svg
export show_lattice_edge_weight_svg
export show_unitcell_svg
export show_diff_lattice_svg
include("lattice/show_lattice.jl")


include("lattice/rspace_matrix.jl")
include("lattice/kmesh.jl")
include("lattice/kspace_matrix.jl")
include("lattice/KP.jl")

export EigenModes, solve_eigenmodes
export DoS
include("lattice/EigenModes.jl")

export qMAT, qMATS
include("lattice/qMAT.jl")


export print
include("lattice/show_datastructures.jl")

#: ------------------------

export ForceConstant
#include("dynmat/bond_bending.jl")
include("dynmat/interatomic_force_matrix.jl")
include("dynmat/ForceConstant.jl")

export to_mat_3x3
export PhonopyCopy
export phonony_from_poscar
export PhonopyCopy_from_vasp
export extract_force_constants

include("dynmat/Phonopy.jl")

export DynamicalMatrixR
export generate_blocks_for_DR
export rspace_dynamical_matrix
export rspace_dynamical_matrix_monolayer
export fold_DymamicalMatrixR
include("dynmat/DynamicalMatrixR.jl")

export DynamicalMatrixQ
export DynamicalMatrix
export ∂DynamicalMatrix∂qα
export consistency_check
export kspace_dynamical_matrix
export kspace_dynamical_matrix_monolayer
export kspace_dynamical_matrix_cylinder
export kspace_dynamical_matrix_ribbon_2D
export dDQda
include("dynmat/DynamicalMatrixQ.jl")

export DynamicalMatricesQ
include("dynmat/DynamicalMatricesQ.jl")

export DynamicalMatricesR
include("dynmat/DynamicalMatricesR.jl")

export DynamicalMatrixSymbolicQ
export kspace_dynamical_matrix_symbolic
export dDQda
include("dynmat/DynamicalMatrixSymbolicQ.jl")

export DynamicalMatrixKp
export DynamicalMatrixSymbolicKp
include("dynmat/DynamicalMatrixKp.jl")


#: ------------------------


export HoppingParameter
export zero_HoppingParameter
export one_onsite_potential, zero_onsite_potential
include("util/MathematicaCForm.jl")
include("hopping/SlaterKosterDict.jl")
include("hopping/HoppingParameter.jl")


export HoppingHamiltonianR
export generate_blocks_for_HR
export rspace_hopping_hamiltonian
export rspace_hopping_hamiltonian_monolayer
export rspace_zero_hopping_hamiltonian
export rspace_zero_hopping_hamiltonian_monolayer
export rspace_onsite_potential
export rspace_onsite_potential_monolayer
export rspace_chemical_potential
export rspace_chemical_potential_monolayer
export solve_eigenmodes!
include("hopping/HoppingHamiltonianR.jl")


export HoppingHamiltonianenR
include("hopping/HoppingHamiltonianenR.jl")

export HoppingHamiltonianQ, HoppingHamiltonian
export ChemicalPotentialQ
export generate_blocks_for_HQ
export kspace_hopping_hamiltonian
export kspace_hopping_hamiltonian_monolayer
export kspace_hopping_hamiltonian_cylinder
export kspace_zero_hopping_hamiltonian
export kspace_zero_hopping_hamiltonian_monolayer
export kspace_zero_hopping_hamiltonian_cylinder
export kspace_chemical_potential
export kspace_chemical_potential_monolayer
export kspace_chemical_potential_cylinder
export dHQda
include("hopping/HoppingHamiltonianQ.jl")

export is_HoppingHamiltonianenQ
export HoppingHamiltonianenQ
include("hopping/HoppingHamiltonianenQ.jl")


export HoppingHamiltonianSymbolicQ
export kspace_hopping_hamiltonian_symbolic
export dHQda
include("hopping/HoppingHamiltonianSymbolicQ.jl")


include("hopping/HoppingHamiltonianKp.jl")
include("hopping/HoppingHamiltonianenKp.jl")


export SuperconductingGapFunctionQ
export SuperconductingGapFunctionenQ
export SuperconductingGapFunction
export kspace_s_wave_Δ
include("hopping/SuperconductingGapFunctionQ.jl")

export BdGHamiltonianQ
export BdGHamiltonianenQ
export BdGHamiltonian
export kspace_BdG_hamiltonian
include("hopping/BdGHamiltonianQ.jl")

export stacking
export ⊕
include("hopping/stacking.jl")

HoppingHamUnion = Union{HoppingHamiltonianQ,HoppingHamiltonianenQ,HoppingHamiltonianKp,HoppingHamiltonianenKp}

#: ------------------------

export generate_kmesh, EnS
include("band/EnS.jl")

export BandStructure
export extract_kpath_to_line
export KpointVecMTuple, VectorKpointVecMTuple, KpointMatMTuple, VectorKpointMatMTuple, KpointVec
include("band/BandStructure.jl")

export BandStructureGrid
include("band/BandStructureGrid.jl")


include("band/band_structure_base.jl")

export band_structure
include("band/bands_no_markers.jl")

export band_structure_with_markers
export band_structure_with_eigenvectors
include("band/bands_with_markers.jl")

export diff_eig, diff_bands
include("band/diff_bands.jl")


export AcousticWeight, OpticalWeight
export LongitudinalProjector, TransverseProjector
export XYProjector, ZProjector
export SublatticeProjector, AtomProjector
export CM_coords
include("band/phonon_projectors.jl")


export compute_band_markers
export compute_band_markers_and_collapse, CBMC
export compute_band_energy_markers_and_collapse, CBEMC
export transform_mat, transform_mat_group
include("band/compute_band_markers.jl")



# include("band/bands_BdG.jl")
# include("band/band_invariants.jl")

export tensor_V, coherent_part_of_kappa, peierls_part_of_kappa
include("thermal/coherent_part_of_kappa.jl")

#: ------------------------


end
