#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse

ROOT = Path(__file__).resolve().parent

sys.path.insert(0, str(ROOT / "mate_module"))
sys.path.insert(0, str(ROOT / "wpnet_module"))

from mate_module.mate.cli.mate import get_opts as get_mate_opts
from mate_module.mate.cli.mate import main as mate_entry
from wpnet_module.gpnet.cli.gpnet import get_opts as get_wpnet_opts
from wpnet_module.gpnet.cli.gpnet import main as wpnet_entry

VERSION = "1.1.0"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    subparsers = parser.add_subparsers()

    parser_mate = subparsers.add_parser("mate", help="MATE module of GPNet")
    get_mate_opts(parser_mate)
    parser_mate.set_defaults(func=mate_entry)

    parser_wpnet = subparsers.add_parser("wpnet", help="WPNET module of GPNet")
    get_wpnet_opts(parser_wpnet)
    parser_wpnet.set_defaults(func=wpnet_entry)

    try:
        args = parser.parse_args()
        args.func(args)
    except AttributeError:
        parser.print_help()


if __name__ == "__main__":
    main()
