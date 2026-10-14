"""Verify the tutorial's three circular platelet families in generated frame 0."""
import argparse
from pathlib import Path
import re
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
folder = parser.parse_args().output_dir


def read_value(path, name):
    text = re.sub(r"#.*", "", path.read_text())
    match = re.search(r"(?:^|;)\s*" + re.escape(name) + r"\s*=([^;]*);", text)
    if match is None:
        raise ValueError(f"Missing {name} in {path}")
    return np.fromstring(match.group(1), sep=" ")


inputs = folder / "inputFiles"
C2G = read_value(inputs / "polycrystal.txt", "C2G1").reshape(3, 3)
F = read_value(inputs / "polycrystal.txt", "F").reshape(3, 3)
b = read_value(inputs / "W.txt", "b_SI")[0]
volume_SI = abs(np.linalg.det(F)) * b**3
lines = (folder / "evl/evl_0.txt").read_text().splitlines()
counts = [int(line) for line in lines[:11]]
start = 11 + sum(counts[:5])
records = np.array([np.fromstring(line, sep=" ") for line in lines[start:start + counts[5]]])
np.testing.assert_allclose(C2G, np.eye(3), atol=1e-12)
offset = 0
for family, normal in [('100', [1, 0, 0]), ('010', [0, 1, 0]), ('001', [0, 0, 1])]:
    inclusion = inputs / f"ellipsoidalInclusionsDensity_{family}.txt"
    target = read_value(inclusion, "targetDensity")[0]
    count = int(np.ceil(target * volume_SI))
    crystal_axes = read_value(inclusion, "crystalAxes").reshape(3, 3)
    means = read_value(inclusion, "axesLognormalDistributionMean_SI")
    crystal_distortion = read_value(inclusion, "crystalEigenDistortion").reshape(3, 3)
    np.testing.assert_allclose(read_value(inclusion, "patternVector_SI"), np.zeros(3))
    np.testing.assert_allclose(read_value(inclusion, "axesLognormalDistributionStd_SI"), np.zeros(3))
    family_records = records[offset:offset + count]
    if len(family_records) != count:
        raise ValueError(f"Missing platelets for family {family}")
    for record in family_records:
        axes = record[4:13].reshape(3, 3) * b
        lengths = np.linalg.norm(axes, axis=0)
        np.testing.assert_allclose(2 * lengths, [40e-9, 40e-9, 3e-9], rtol=1e-12, atol=1e-22)
        np.testing.assert_allclose(axes[:, 2] / lengths[2], normal, atol=1e-12)
        np.testing.assert_allclose(axes, C2G @ crystal_axes @ np.diag(means), rtol=1e-12, atol=1e-22)
        np.testing.assert_allclose(record[13:22].reshape(3, 3),
                                   C2G @ crystal_distortion @ C2G.T, atol=1e-14)
    offset += count
    print(f"PASS: {count} circular platelets on ({family}), diameter 40 nm, thickness 3 nm.")
if counts[5] != offset:
    raise ValueError(f"Expected {offset} total platelets, found {counts[5]}")
print("PASS: all three families have correct geometry/tensors and no positional pattern.")
