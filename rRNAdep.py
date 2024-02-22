import os
import sys

from libs import cli


if __name__ == '__main__':
    args = cli.parse_args()
    args.func(args)
