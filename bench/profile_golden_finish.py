"""Profile one checkpoint-completion attempt; retain data on interruption."""
import argparse
import cProfile
from pathlib import Path
from complete_golden_campaign import finish

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--run',type=Path,required=True)
parser.add_argument('--index',type=int,required=True)
parser.add_argument('--cut-workers',type=int,default=1)
args=parser.parse_args()
profile=cProfile.Profile();profile.enable()
try:
    finish(args)
finally:
    profile.disable()
    profile.dump_stats(str(args.run/str(args.index)/'completion_profile.pstats'))
