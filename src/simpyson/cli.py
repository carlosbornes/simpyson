import argparse
from simpyson.gui import main as gui_main

def main():
    parser = argparse.ArgumentParser(prog="simpyson")
    subparsers = parser.add_subparsers(dest="command")

    parser_gui = subparsers.add_parser("gui", help="Launch the GUI")
    parser_gui.set_defaults(func=lambda args: gui_main())

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

