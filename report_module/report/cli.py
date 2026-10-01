def get_opts(parser):
    parser.add_argument("-i", "--input", help="Input mat or select file", required=True)
    parser.add_argument("-m", "--mat", help="Input population mat file", required=True)
    parser.add_argument("-o", "--output", help="Output report directory", required=True)
