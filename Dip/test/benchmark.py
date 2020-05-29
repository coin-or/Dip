import unittest
from ipet import Experiment
from ipet.evaluation import IPETEvaluation
import os
import logging
import subprocess
import time
import git


logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

repo = git.Repo(search_parent_directories=True)
sha = repo.head.object.hexsha

DIP_EXE = os.environ.get("DIP_EXE", None)
if not DIP_EXE:
    raise Exception(
        "Env variable DIP_EXE not set. Do 'DIP_EXE={path_to_exe} python ...'"
    )
if not os.path.isfile(DIP_EXE):
    raise Exception(f"Env variable DIP_EXE file does not exist {DIP_EXE}")

DATA_DIR = os.environ.get("DATA_DIR", None)
if not DATA_DIR:
    raise Exception(
        "Env variable DATA_DIR not set. Do 'DATA_DIR={path_to_data} python ...'"
    )
if not os.path.isdir(DATA_DIR):
    raise Exception(f"Env variable DATA_DIR file does not exist {DATA_DIR}")

srcdir = os.path.dirname(os.path.realpath(__file__))


class TestBenchmark(unittest.TestCase):
    def do_test(self, instancefile, blockfile, blockformat):
        outfile = os.path.splitext(os.path.basename(instancefile))[0] + ".out"
        solufile = os.path.join(srcdir, "benchmark.solu")

        # need something here to make it work wit makefiles
        output = subprocess.Popen(
            [
                DIP_EXE,
                "--BlockFileFormat",
                blockformat,
                "--Instance",
                instancefile,
                "--BlockFile",
                blockfile,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        stdout, stderr = output.communicate()

        # dump output
        with open(outfile, "w") as f:
            f.write(f"@01 {instancefile}\n")
            f.write(f"GitHash: {sha}\n")
            f.write(f"@03 {time.time()}\n")
            f.write(stdout.decode("utf-8"))
            f.write("\n")
            f.write(f"@04 {time.time()}\n")
            f.write(f"=ready=\n")

        experiment = Experiment()
        experiment.addOutputFile(outfile)
        # experiment.addSoluFile(solufile)
        experiment.collectData()

        tr = experiment.getTestRuns()[0]
        status = tr.data.iloc[0]["Status"]

        # maybe assert that the running time has not gone through the roof by some mistake
        # throuh ipet? Or just using a local file, will depend on machine, settings, etc.
        solvingTime = tr.data.iloc[0]["SolvingTime"]

        return status, solvingTime

    def test_retail3(self):
        instancefile = os.path.join(DATA_DIR, "retail3.mps")
        blockfile = os.path.join(DATA_DIR, "retail3.block")
        blockformat = "Pair"

        status, solvingTime = self.do_test(instancefile, blockfile, blockformat)

        assert status == "ok"

    def test_atm_5_10_1(self):
        instancefile = os.path.join(DATA_DIR, "atm_5_10_1.mps")
        blockfile = os.path.join(DATA_DIR, "atm_5_10_1.block")
        blockformat = "List"

        status, solvingTime = self.do_test(instancefile, blockfile, blockformat)

        assert status == "ok"

    def test_wedding_16(self):
        instancefile = os.path.join(DATA_DIR, "wedding_16.mps")
        blockfile = os.path.join(DATA_DIR, "wedding_16.block")
        blockformat = "Pair"

        status, solvingTime = self.do_test(instancefile, blockfile, blockformat)

        assert status == "ok"

    def test_block_milp(self):
        instancefile = os.path.join(DATA_DIR, "block_milp.lp")
        blockfile = os.path.join(DATA_DIR, "block_milp.dec")
        blockformat = "ZIBList"

        status, solvingTime = self.do_test(instancefile, blockfile, blockformat)

        assert status == "ok"


if __name__ == "__main__":
    unittest.main()
