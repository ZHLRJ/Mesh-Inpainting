# Neural MeshImplant: Data-Free Learning for Mesh Repair

Official implementation of the paper:

**Neural MeshImplant: Data-Free Learning for Mesh Repair**

---

## Overview
![teaser image](media/Teaser.png)
This repository implements a **data-free neural mesh inpainting framework** that restores missing regions directly on the input mesh while strictly preserving its original topology and connectivity.

Unlike traditional dataset-driven methods, our approach:

- Operates on a **single incomplete mesh**
- Requires **no external training datasets**
- Preserves original vertex indexing and mesh structure
- Supports direct use in downstream applications

The framework consists of two stages:

1. **Topology-Preserving Geometric Implantation**
2. **Self-Supervised Neural Refinement**

---

## Key Contributions

- ✔ Data-free mesh repair
- ✔ Geodesic-based cut-and-stitch hole filling
- ✔ Topology-preserving surface implantation
- ✔ Self-supervised GCN refinement
- ✔ Compatible with PDE solvers and FEM pipelines
- ✔ Scalable to high-resolution meshes

---

## Method Overview

### 1️⃣ Geometric Stage (Implantation)

We first preprocess the incomplete mesh to construct a watertight auxiliary surface using:

- Zero-level isosurface extraction (Signed Heat Method)
- Poisson surface reconstruction
- QSLIM simplification
- Triangular remeshing

Hole boundaries are projected onto this surface, and smooth **geodesic cut paths** are generated.  
We then stitch patches while preserving manifoldness using a half-edge structure.

This stage guarantees:

- Topological consistency
- Edge-manifoldness
- Parallel processing of multiple holes

---

### 2️⃣ Self-Supervised Neural Refinement

We adopt a mesh autoencoder architecture with:

- QSLIM-based pooling / unpooling
- Graph Convolutional Networks
- Vertex displacement prediction
- Laplacian smoothing initialization

Loss functions:

- Vertex position loss
- Face normal consistency loss
- Bilateral normal filtering regularization

This stage fine-tunes repaired regions without requiring any external dataset.

---

## Applications

Repaired meshes can be directly used for:

- ✔ Surface PDE continuation
- ✔ Logarithmic map computation
- ✔ Parallel transport
- ✔ FEM simulation
- ✔ 3D printing
- ✔ Point cloud completion

Because connectivity and vertex indexing are preserved,  
downstream solvers can reuse boundary conditions without remeshing.

---

## Installation

Clone repository:

```bash
git clone --recursive https://github.com/YOUR_USERNAME/Mesh-Inpainting.git
