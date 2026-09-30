#!/usr/bin/env python3
from mate.cli.mate import get_opts
from mate.cli.mate import main
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    get_opts(parser)
    opts = parser.parse_args()
    main(opts)
