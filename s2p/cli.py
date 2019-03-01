import os
import shutil
import argparse

import s2p


def main():
    """
    Command line parsing for s2p command line interface.
    """
    parser = argparse.ArgumentParser(description=('S2P: Satellite Stereo '
                                                  'Pipeline'))
    parser.add_argument('config', metavar='config.json',
                        help=('path to a json file containing the paths to '
                              'input and output files and the algorithm '
                              'parameters'))
    args = parser.parse_args()

    user_cfg = s2p.read_config_file(args.config)

    s2p.main(user_cfg)

    # Backup input file for sanity check
    if not args.config.startswith(os.path.abspath(s2p.cfg['out_dir']+os.sep)):
        shutil.copy2(args.config,os.path.join(s2p.cfg['out_dir'],'config.json.orig'))
