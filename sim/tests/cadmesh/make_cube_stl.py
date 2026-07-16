#!/usr/bin/env python3
"""Write a watertight ASCII-STL cube (outward CCW normals) for the CADMesh gate.

CADMesh's built-in reader is ASCII-STL only (binary needs ASSIMP, which we do
not link), so the S4 geometry generator must emit ASCII STL. This is the minimal
shared writer used both by the dependency test and as the reference format.

    python make_cube_stl.py cube_ascii.stl [side_mm=10]
"""
import sys


def write_cube_ascii(path: str, side: float = 10.0) -> float:
    h = side / 2.0
    v = [(-h, -h, -h), (h, -h, -h), (h, h, -h), (-h, h, -h),
         (-h, -h, h), (h, -h, h), (h, h, h), (-h, h, h)]
    quads = [(0, 3, 2, 1, (0, 0, -1)), (4, 5, 6, 7, (0, 0, 1)),
             (0, 1, 5, 4, (0, -1, 0)), (2, 3, 7, 6, (0, 1, 0)),
             (1, 2, 6, 5, (1, 0, 0)), (0, 4, 7, 3, (-1, 0, 0))]
    tris = []
    for a, b, c, d, n in quads:
        tris += [(n, v[a], v[b], v[c]), (n, v[a], v[c], v[d])]
    with open(path, "w") as f:
        f.write("solid cube\n")
        for n, p1, p2, p3 in tris:
            f.write(f"  facet normal {n[0]:.1f} {n[1]:.1f} {n[2]:.1f}\n")
            f.write("    outer loop\n")
            for p in (p1, p2, p3):
                f.write(f"      vertex {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            f.write("    endloop\n  endfacet\n")
        f.write("endsolid cube\n")
    return side ** 3


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "cube_ascii.stl"
    side = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
    vol = write_cube_ascii(out, side)
    print(f"wrote {out}: side={side} mm, expected volume={vol:.0f} mm^3")
