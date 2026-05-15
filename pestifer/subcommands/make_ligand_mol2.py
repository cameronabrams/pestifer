# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
The make-ligand-mol2 subcommand: generate CGenFF-ready mol2 files for
every HETATM ligand in a PDB whose resname is not already defined in the
loaded CHARMM force field.

The resulting mol2 files are intended as input either to a local
``cgenff`` binary or to the CGenFF web tool (``cgenff.com``). The same
internals are reused by the build workflow when it encounters an
unknown ligand.
"""
from __future__ import annotations

import argparse as ap
import logging
from dataclasses import dataclass
from pathlib import Path

from . import Subcommand
from ..charmmff.ligand_paramgen.orchestrator import (
    LigandGenSummary,
    fetch_or_load_pdb,
    generate_ligand_mol2s,
)
from ..core.config import Config

logger = logging.getLogger(__name__)


@dataclass
class MakeLigandMol2Subcommand(Subcommand):
    name: str = "make-ligand-mol2"
    short_help: str = (
        "generate CGenFF-ready mol2 files for unknown HETATM ligands in a PDB"
    )
    long_help: str = (
        "For each HETATM residue in the input PDB whose resname is not "
        "covered by the loaded CHARMM force-field topologies, fetch a "
        "SMILES from the RCSB Chemical Component Dictionary, protonate "
        "the ligand at the chosen pH (default 7.4) using RDKit + "
        "Dimorphite-DL, and write a Tripos mol2 ready to feed the CGenFF "
        "binary or web tool. SOURCE may be a 4-letter RCSB PDB code or a "
        "path to a local .pdb file."
    )
    func_returns_type: type = dict

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        overrides: dict[str, str] = {}
        for kv in args.smiles or []:
            if "=" not in kv:
                raise ValueError(
                    f"--smiles expects RESNAME=SMILES, got {kv!r}"
                )
            resname, smi = kv.split("=", 1)
            overrides[resname.strip().upper()] = smi.strip()

        config = Config().configure_new()
        cc = config.RM.charmmff_content

        parsed = fetch_or_load_pdb(args.source)
        summary: LigandGenSummary = generate_ligand_mol2s(
            parsed,
            cc,
            outdir=Path(args.outdir),
            ph=args.ph,
            smiles_overrides=overrides,
        )

        _print_summary(summary, args.outdir)
        return {
            "successes": [r.resname for r in summary.successes],
            "failures": {r.resname: r.status for r in summary.failures},
            "skipped_known": summary.skipped_known,
        }

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument(
            "source",
            type=str,
            help="4-letter RCSB PDB code, or path to a local .pdb file",
        )
        self.parser.add_argument(
            "--outdir",
            type=str,
            default="ligand_mol2",
            help="directory to write {RESNAME}.mol2 files into "
                 "(default: %(default)s)",
        )
        self.parser.add_argument(
            "--ph",
            type=float,
            default=7.4,
            help="target pH for ligand protonation (default: %(default)s)",
        )
        self.parser.add_argument(
            "--smiles",
            action="append",
            metavar="RESNAME=SMILES",
            help="override the SMILES used for a given resname; repeatable. "
                 "Useful when RCSB has no entry, or to force a specific "
                 "tautomer/charge state.",
        )
        return self.parser


def _print_summary(summary: LigandGenSummary, outdir: str) -> None:
    n_ok = len(summary.successes)
    n_fail = len(summary.failures)
    print()
    print(f"=== make-ligand-mol2 summary (outdir: {outdir}) ===")
    if summary.skipped_known:
        print(
            f"Already in CHARMM ({len(summary.skipped_known)}): "
            f"{', '.join(summary.skipped_known)}"
        )
    if n_ok:
        print(f"\nGenerated {n_ok} mol2 file(s):")
        for r in summary.successes:
            print(f"  {r.resname:>5s}  ->  {r.path}")
    if n_fail:
        print(f"\nFailed ({n_fail}):")
        for r in summary.failures:
            print(f"  {r.resname:>5s}  [{r.status}]  {r.message}")
    if not n_ok and not n_fail:
        print("No unknown ligands found.")
