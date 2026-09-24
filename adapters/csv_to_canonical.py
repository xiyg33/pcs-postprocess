"""Example: convert a partner CSV with canonical column names to MAT/JSON."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import savemat


def convert(csv_path: Path, meta_path: Path, output_dir: Path) -> Path:
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    with csv_path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        names = reader.fieldnames
        if not names:
            raise ValueError("CSV needs a header row")
        data = np.array([[float(row[name]) for name in names] for row in reader], dtype=float)
    if data.ndim != 2 or data.shape[0] < 2:
        raise ValueError("CSV needs at least two data rows")
    header = np.empty((len(names), 2), dtype=object)
    for index, name in enumerate(names):
        header[index] = (name, "")
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"case_{int(meta['case_index']):04d}_{meta['case_name']}"
    output = output_dir / f"{stem}.mat"
    savemat(output, {"data": data, "header": header}, do_compression=True)
    (output_dir / f"{stem}_meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--meta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(convert(args.csv, args.meta, args.output))
