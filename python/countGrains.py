"""Count mesh regions using the same element tags as MoDELib's SimplexReader.

Usage (from the repository root):
    python3 python/countGrains.py Library/Meshes/poly30_periodicXYZ.msh

Supports ASCII Gmsh 2.x files with linear or quadratic simplices. The highest
simplex dimension present determines whether regions are volumes, areas, or
lines. Boundary-element tags are excluded from the grain count.
"""

import argparse
from pathlib import Path


def countGrains(meshFileName):
    """Return the number of distinct region IDs in a MoDELib mesh file."""
    dimensions = {1: 1, 8: 1, 2: 2, 9: 2, 4: 3, 11: 3}
    nodeCounts = {1: 2, 8: 3, 2: 3, 9: 6, 4: 4, 11: 10}
    regions = {1: set(), 2: set(), 3: set()}
    haveFormat = False
    haveElements = False
    with Path(meshFileName).open('r', encoding='ascii') as mesh:
        for line in mesh:
            if line.strip() == '$MeshFormat':
                fields = mesh.readline().split()
                if len(fields) != 3 or not fields[0].startswith('2.') or fields[1] != '0':
                    raise ValueError('Expected an ASCII Gmsh 2.x mesh')
                if mesh.readline().strip() != '$EndMeshFormat':
                    raise ValueError('Invalid MeshFormat section')
                haveFormat = True
            elif line.strip() == '$Elements':
                if not haveFormat:
                    raise ValueError('Missing MeshFormat before Elements')
                count = int(mesh.readline())
                if count < 0:
                    raise ValueError('Negative element count')
                for _ in range(count):
                    fields = [int(value) for value in mesh.readline().split()]
                    if len(fields) < 3:
                        raise ValueError('Incomplete element record')
                    elementID, elementType, tagCount = fields[:3]
                    if tagCount < 0 or len(fields) < 3 + tagCount:
                        raise ValueError(f'Invalid tags for element {elementID}')
                    if elementType in dimensions:
                        if tagCount < 2 or len(fields) != 3 + tagCount + nodeCounts[elementType]:
                            raise ValueError(f'Element {elementID} is missing vertices or region tag')
                        regions[dimensions[elementType]].add(fields[4])
                if mesh.readline().strip() != '$EndElements':
                    raise ValueError('Invalid Elements section')
                haveElements = True
    if not haveElements:
        raise ValueError('Missing Elements section')
    for dimension in (3, 2, 1):
        if regions[dimension]:
            return len(regions[dimension])
    raise ValueError('No supported simplex elements found')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('meshFileName', help='Path to an ASCII Gmsh 2.x mesh')
    args = parser.parse_args()
    try:
        print(countGrains(args.meshFileName))
    except (OSError, ValueError) as error:
        parser.exit(1, f'{parser.prog}: {error}\n')


if __name__ == '__main__':
    main()
