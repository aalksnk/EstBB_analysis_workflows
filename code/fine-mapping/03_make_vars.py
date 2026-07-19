#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np

# compatibility for packages still using np.int / np.float
if not hasattr(np, "int"):
    np.int = int
if not hasattr(np, "float"):
    np.float = float 
from ldstore.bcor import bcor  # from the ldstore package

def main():
    parser = argparse.ArgumentParser(
        description="Extract SNP ID order from .bcor files into .vars files."
    )
    parser.add_argument(
        "--bcor-dir",
        required=True,
        help="Directory containing .bcor files."
    )
    parser.add_argument(
        "--id-column",
        default="rsid",
        help="Column in bcor metadata to use as variant ID (default: rsid)."
    )
    args = parser.parse_args()

    bcor_dir = Path(args.bcor_dir)
    if not bcor_dir.is_dir():
        raise SystemExit(f"{bcor_dir} is not a directory")

    bcor_files = sorted(bcor_dir.glob("*.bcor"))
    if not bcor_files:
        raise SystemExit(f"No .bcor files found in {bcor_dir}")

    for bcor_path in bcor_files:
        print(f"=== {bcor_path.name} ===")
        bc = bcor(str(bcor_path))

        meta = bc.getMeta()
        if args.id_column not in meta.columns:
            raise SystemExit(
                f"{bcor_path}: metadata missing '{args.id_column}' column. "
                f"Available: {list(meta.columns)}"
            )

        ids = meta[args.id_column].astype(str).tolist()
        out_path = bcor_path.with_suffix(".vars")

        with out_path.open("w") as out:
            out.write("\n".join(ids) + "\n")

        print(f"  Wrote {len(ids)} IDs -> {out_path}")

if __name__ == "__main__":
    main()
