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

#Define parameters of the RBF sampler
rSampler = 0.8 

# Setup time steps
steps = int(2)

# Setup output interval
outputInterval  = int(1)

# Spatial/temporal samplers
t_sampler = mui4py.TemporalSamplerExact()
s_sampler = mui4py.SamplerPseudoNearestNeighbor(rSampler)

# Commit ZERO step
MUI_Interfaces["MLDataInterface"].barrier(0)

# Fetch data from the other solver
peer_size_fetch = MUI_Interfaces["MLDataInterface"].fetch(name_fetch_parameter)

print("{Python} peer_size: ", peer_size, " peer_size_fetch: ", peer_size_fetch, "\n")

assert peer_size == peer_size_fetch

name_list_fetch_cellSize = [f"cellSize_{r}" for r in range(peer_size_fetch)]

# Loop over all peer ranks and generate their names
for r in range(peer_size_fetch):
    data_types_fetch_cellSize = {name_list_fetch_cellSize[r]: mui4py.INT64}
    MUI_Interfaces["MLDataInterface"].set_data_types(data_types_fetch_cellSize)
    print(f"Rank {r}: prepared name for rank {r}: {name_list_fetch_cellSize[r]}")

# Fetch data from the other solver
peer_size_fetch = [MUI_Interfaces["MLDataInterface"].fetch(name) for name in name_list_fetch_cellSize]

print(f"peer_size_fetch: = {peer_size_fetch}")

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

    #assert Npoints > 0
