"""
/*****************************************************************************
* Parallel Partitioned Multi-Physics Simulation Framework (ParaSiF)          *
*                                                                            *
* Copyright (C) 2025 The ParaSiF Development Team                            *
* All rights reserved                                                        *
*                                                                            *
* This software is licensed under the GNU General Public License version 3   *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/
#
# @file dummyPython.py
# @author W. Liu
# @date 09 October 2025
# @brief Two-way date exchange between OpenFOAM and a dummy python solver using MUI library
#
"""

from mpi4py import MPI
import mpi4py
import datetime
import numpy as np
import time
import os

# Include MUI header file and configure file 
import mui4py

# Common world claims 
MUI_COMM_WORLD = mui4py.mpi_split_by_app()

# Declare MPI ranks
rank = MUI_COMM_WORLD.Get_rank()

# Declare MPI size
size = MUI_COMM_WORLD.Get_size()

peer_size = MPI.COMM_WORLD.Get_size() - size

# Define MUI dimension
dimensionMUI = 3

# Define the name of push/fetch values
name_push = "prghDouble"
name_fetch  = "prgh"
name_fetch_parameter  = "nProcs"

# Define MUI push/fetch data types
data_types_push = {name_push: mui4py.FLOAT64}
data_types_fetch = {name_fetch: mui4py.FLOAT64,
				    name_fetch_parameter: mui4py.INT64} 

# MUI interface creation
domain = "python"
config3d = mui4py.Config(dimensionMUI, mui4py.FLOAT64)
iface = ["MLDataInterface"]
MUI_Interfaces = mui4py.create_unifaces(domain, iface, config3d) 
MUI_Interfaces["MLDataInterface"].set_data_types(data_types_fetch)
MUI_Interfaces["MLDataInterface"].set_data_types(data_types_push)

print("mpi4py.get_config(): ", mpi4py.get_config(), "\n")

# Define the forget steps of MUI to reduce the memory
forgetSteps = int(5)
# Define the synchronised boolen for MUI smart send
synchronised=False

#Define parameters of the tol sampler
tol = 1e-3

# Setup time steps
time_step = int(2)
nOuterCorrectors = int(3)
nCorrectors = int(4)
nNonOrthogonalCorrectors = int(6)
steps = int(time_step * nOuterCorrectors * nCorrectors * nNonOrthogonalCorrectors)

# Setup output interval
outputInterval  = int(1)

# Spatial/temporal samplers
t_sampler = mui4py.TemporalSamplerExact()
s_sampler = mui4py.SamplerExact(tol)

# Barrier zero step
MUI_Interfaces["MLDataInterface"].barrier(0)

# Fetch data from the other solver
peer_size_fetch = MUI_Interfaces["MLDataInterface"].fetch(name_fetch_parameter)

print("{Python} peer_size: ", peer_size, " peer_size_fetch: ", peer_size_fetch, "\n")

assert peer_size == peer_size_fetch

name_list_fetch_cellSize = [f"cellSize_{r}" for r in range(peer_size_fetch)]
name_list_fetch_CoordX = [f"cellX_{r}" for r in range(peer_size_fetch)]
name_list_fetch_CoordY = [f"cellY_{r}" for r in range(peer_size_fetch)]
name_list_fetch_CoordZ = [f"cellZ_{r}" for r in range(peer_size_fetch)]

# Loop over all peer ranks and generate their names
for r in range(peer_size_fetch):
    data_types_fetch_cellSize = {name_list_fetch_cellSize[r]: mui4py.INT64}
    MUI_Interfaces["MLDataInterface"].set_data_types(data_types_fetch_cellSize)
    data_types_fetch_CoordX = {name_list_fetch_CoordX[r]: mui4py.FLOAT64}
    MUI_Interfaces["MLDataInterface"].set_data_types(data_types_fetch_CoordX)
    data_types_fetch_CoordY = {name_list_fetch_CoordY[r]: mui4py.FLOAT64}
    MUI_Interfaces["MLDataInterface"].set_data_types(data_types_fetch_CoordY)
    data_types_fetch_CoordZ = {name_list_fetch_CoordZ[r]: mui4py.FLOAT64}
    MUI_Interfaces["MLDataInterface"].set_data_types(data_types_fetch_CoordZ)
    print(f"Rank {r}: prepared name for rank {r}: {name_list_fetch_cellSize[r]}")

# Fetch data from the other solver
cell_size_fetch = [MUI_Interfaces["MLDataInterface"].fetch(name) for name in name_list_fetch_cellSize]

print(f"cell_size_fetch: = {cell_size_fetch}")

coords_all = np.empty(peer_size_fetch, dtype=object)

for r in range(peer_size_fetch):

    locp = np.column_stack((np.arange(0, int(cell_size_fetch[r])), np.zeros(cell_size_fetch[r]), np.zeros(cell_size_fetch[r])))
    print("{Python} locp: ", locp)
    N = int(cell_size_fetch[r])  # ensure it’s an integer count

    coords_all[r] = np.zeros((N, 3))

    coords_x = MUI_Interfaces["MLDataInterface"].\
                            fetch_many(name_list_fetch_CoordX[r],
                            locp,
                            0,
                            s_sampler,
                            t_sampler)
    coords_y = MUI_Interfaces["MLDataInterface"].\
                            fetch_many(name_list_fetch_CoordY[r],
                            locp,
                            0,
                            s_sampler,
                            t_sampler)
    coords_z = MUI_Interfaces["MLDataInterface"].\
                            fetch_many(name_list_fetch_CoordZ[r],
                            locp,
                            0,
                            s_sampler,
                            t_sampler)

    coords_all[r] = np.column_stack((coords_x, coords_y, coords_z))

    print(f"Rank {r}: prepared name for rank {r}: {name_list_fetch_CoordX[r]}, {name_list_fetch_CoordY[r]}, {name_list_fetch_CoordZ[r]}")

print([coords.shape for coords in coords_all])
print(coords_all)

# Output the initial pseudo scalar field
output_dir = "coupling_results"
if rank == 0:
    os.makedirs(output_dir, exist_ok=True)

outputFileName = os.path.join(output_dir, "scalar_field.csv")
with open(outputFileName, "w+") as outputFile:
    outputFile.write("\"X\",\"Y\",\"Z\",\"scalar_field\"\n")
outputFile.close() 

# Begin time loops
for t in range(1, (steps+1)):

    print("\n")
    print("{Python} ", t," Step \n")

    prgh_all = np.empty(peer_size_fetch, dtype=object)

    for r in range(peer_size_fetch):
        prgh_all[r] = MUI_Interfaces["MLDataInterface"].\
                            fetch_many(name_fetch,
                            coords_all[r],
                            t,
                            s_sampler,
                            t_sampler)
    print(prgh_all)
    print(type(prgh_all))
    print("dtype:", prgh_all.dtype)
    print("shape:", prgh_all.shape)
    print("ndim:", prgh_all.ndim)

    # Push data to the other solver
    prghDouble = np.array([2.0 * arr for arr in prgh_all], dtype=object)
    for r in range(peer_size_fetch):
        MUI_Interfaces["MLDataInterface"].push_many(name_push, coords_all[r], prghDouble[r])

    # Commit 't' step of MUI
    MUI_Interfaces["MLDataInterface"].commit(t)
