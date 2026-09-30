#!/usr/bin/env python3
from gpnet.cli.gpnet import main
from gpnet.cli.gpnet import get_opts
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    get_opts(parser)
    opts = parser.parse_args()
    main(opts)
