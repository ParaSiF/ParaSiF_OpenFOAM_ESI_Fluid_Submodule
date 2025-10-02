# ParaSiF OpenFOAM (ESI) Fluid Solver

This repository contains the **OpenFOAM (ESI) solver** integrated with the ParaSiF Parallel Partitioned Simulation Framework.
It is maintained as a **submodule** of the main ParaSiF repository: [ParaSiF Main Repository](https://github.com/ParaSiF/ParaSiF).

---

## Overview

The OpenFOAM (ESI) Fluid Solver in ParaSiF allows the simulation of *fluid domains* in **multi-physics partitioned simulations**.
It is designed to interface with other solvers (e.g., structural solvers) via the **[MUI coupling library](https://mxui.github.io/)**.

Key features:

- Two-way partitioned coupling with structural solvers.
- Supports parallel execution on HPC systems.
- Modular design: can be replaced or updated without affecting other solvers.
- Compatible with precompiled or user-installed OpenFOAM versions.

---

## Compatible Codebase

This solver has been tested and is compatible with **[OpenFOAM v2506](https://www.openfoam.com/news/main-news/openfoam-v2506)**.

> Users are recommended to use this version to ensure full compatibility with ParaSiF OpenFOAM (ESI) solvers.

---

## Location in the Main ParaSiF Repository

`ParaSiF/src/fluid/OpenFOAM_ESI/`

---

## Repository Structure

```
ParaSiF/src/fluid/OpenFOAM_ESI/
├── doc/                     # Documentation folder
├── src/                     # ParaSiF-specific source code folder
│ ├── applications/          # ParaSiF-specific applications
│ │ └── solvers/             # ParaSiF-specific solvers
│ │  ├── pimpleFSIFoam       # source code of pimpleFSIFoam solver
│ │  └── interFSIFoam        # source code of interFSIFoam solver
│ └── libs/                  # Libs used by solvers
└── test/                    # test folder
 ├── mui_communication_test  # Tests for MUI communication under OpenFOAM codebase
 ├── pimpleFSIFoam_test      # Test cases for pimpleFSIFoam solver
 └── interFSIFoam_test       # Test cases for interFSIFoam solver
```
## Installation

**Note:** This solver is a submodule of ParaSiF. Follow the main ParaSiF repository instructions to initialise submodules and install global dependencies.

### Steps

1. **Obtain and install the codebase**

   - Initialise the submodule from the main ParaSiF repository (or clone this repository).
   - Ensure the correct version of OpenFOAM is installed.

2. **Source your OpenFOAM environment**

```bash
source /path/to/OpenFOAM/etc/bashrc
```

3. **Compile the ParaSiF libraries (if applicable)**

```bash
cd src/libs
wmake libso <lib_name>
```

4. **Compile the ParaSiF solvers**

```bash
cd src/applications/solvers/
wmake <solver_name>
```

> If you encounter errors related to mui.h, update the paths in the Make/options files accordingly and rerun wmake.

## Running Tests and Example Cases

Test cases are located in the test/ folder:

To run a test:

```bash
source /path/to/OpenFOAM/etc/bashrc
cd test/XXX_test
./Allrun
```

Output will be saved in runData/ for analysis.

For integrated example cases with other solvers, see the example/ folder in the main ParaSiF repository.

## Contributing

ParaSiF, including this submodule, is an **open-source project**, and contributions from the community are warmly welcomed.

There are many ways you can help improve this submodule, including:

- Adding new features, libs or solvers
- Improving documentation, tests and examples
- Fixing bugs or refining existing functionality
- Sharing feedback and suggestions for enhancements

Your contributions, whether large or small, are highly valued and help make ParaSiF a stronger resource for the research community.

For detailed guidance on contributing, please see the [CONTRIBUTING.md](https://github.com/ParaSiF/ParaSiF/blob/main/CONTRIBUTING.md) in the main ParaSiF repository.

## License

Copyright (C) 2021–2025 The ParaSiF Development Team.  
Licensed under the **GNU General Public License v3 (GPL-3.0)**.

## Contact

For questions or contributions, please contact the ParaSiF Development Team
