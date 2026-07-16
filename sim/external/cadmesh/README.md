# Vendored CADMesh (v2.0.3)

Single-header STL/PLY/OBJ → `G4TessellatedSolid` loader for the S4 lattice
geometries. Source: https://github.com/christopherpoole/CADMesh (MIT).

**Header-only, no external deps for STL.** Just `#include "CADMesh.hh"` and add
this dir to the include path. The built-in reader is the default
(`CADMESH_DEFAULT_READER = BuiltIn`); ASSIMP and tetgen are behind
`USE_CADMESH_ASSIMP_READER` / `USE_CADMESH_TETGEN` (left undefined → not needed).

Usage:
```cpp
auto mesh = CADMesh::TessellatedMesh::FromSTL("part.stl");
G4TessellatedSolid* solid = mesh->GetTessellatedSolid();  // already SetSolidClosed(true)
```

Verified to build + link + round-trip a cube against Geant4 11.3.2 (S3.5,
see CADMESH_CHECK.md).
