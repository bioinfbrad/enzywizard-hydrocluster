from __future__ import annotations
from argparse import Namespace, ArgumentParser
from ..services.hydrocluster_service import run_hydrocluster_service


def add_hydrocluster_parser(parser: ArgumentParser) -> None:
    parser.add_argument("-i","--input_path", required=True, help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to the output directory for saving the JSON report.")
    parser.add_argument("--cutoff", type=float, default=10.0, help="Minimum residue-residue contact area cutoff for hydrophobic cluster connection.")
    parser.set_defaults(func=run_hydrocluster)

def run_hydrocluster(args: Namespace) -> None:
    run_hydrocluster_service(input_path=args.input_path, output_dir=args.output_dir, cutoff_area=args.cutoff)

