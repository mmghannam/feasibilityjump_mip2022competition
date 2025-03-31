import sys,os
import argparse
import gurobipy as gp

parser = argparse.ArgumentParser(
    description="Gurobi format conversion and presolving."
)
parser.add_argument(
    "inputfile",
    help="Input file (MPS, LP, etc.)."
)
parser.add_argument(
    "outputfile",
    help="Output file (MPS, LP, etc.)."
)
parser.add_argument(
    "--presolve",
    action="store_true",
    help="Apply presolve before saving."
)
parser.add_argument(
    "--no-overwrite",
    action="store_true",
    help="Abort if output file already exists."
)

args = parser.parse_args()

if args.no_overwrite and os.path.exists(args.outputfile):
    print("file exists.")
    sys.exit()


env = gp.Env()
env.setParam("Threads", 1)
#env.setParam("Aggregate", 0)
#env.setParam("Presolve", 2)
#env.setParam("PreSOS1Encoding", 0)
#env.setParam("PreSOS2Encoding", 0)
#env.setParam("PrePasses", 1)

model = gp.read(args.inputfile, env)

if args.presolve:
    model_presolved = model.presolve()
    model_presolved.write(args.outputfile)
else:
    model.write(args.outputfile)
