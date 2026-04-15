from __future__ import annotations

import argparse

from .commands.hydrocluster import add_hydrocluster_parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzywizard-hydrocluster",
        description="EnzyWizard-HydroCluster: Identify hydrophobic clusters from input CIF/PDB file and generate a detailed JSON report."
    )
    add_hydrocluster_parser(parser)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)