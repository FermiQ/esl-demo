#!/usr/bin/env python3


import argparse
import os
import sys
import shlex, subprocess

from time import gmtime, strftime
from ruamel.yaml import YAML


# Process command-line arguments
parser = argparse.ArgumentParser(
        description="Checks ESL Demo test results against references")
parser.add_argument("-a", "--abs-tol", type=float, default=1.0e-7, 
        help="Absolute tolerance for comparisons (overrides configuration)")
parser.add_argument("-c", "--config", default="esl-testsuite.yml",
        help="Config file containing the SIESTA test suite specifications")
parser.add_argument("-d", "--details", action="store_true", default=False,
        help="Display detailed information about each test")
parser.add_argument("-o", "--output", default="tests-report.yml",
        help="File to store a YAML version of the report with full details")
parser.add_argument("-r", "--refdir", default="YAML_Refs",
        help="Directory containing reference YAML files")
parser.add_argument("-t", "--testdir", default=".",
        help="Directory containing the test cases with their outputs")
parser.add_argument("-s", "--single", default="ag.yml",
        help="Single test output yaml")
parser.add_argument("-e", "--exe", default=".",
        help="executable to be tested")
parser.add_argument("-y", "--yaml", default="output.yml",
        help="Relative path where each test stores its YAML data")
args = parser.parse_args()

# Check command-line arguments
if ( not os.path.exists(args.config) ):
    parser.error("config file not found: '%s'" % args.config)
if ( not os.path.isdir(args.refdir) ):
    parser.error("reference directory not found: '%s'" % args.refdir)
if ( not os.path.isdir(args.testdir) ):
    parser.error("tests directory not found: '%s'" % args.testdir)
if ( not os.path.exists(args.exe) ):
    parser.error("executable does not exist: '%s'" % args.exe)

print("# Running {0:s} {1:s}.yml > output.yml".format(args.exe,args.single))
cmd="{0:s} {1:s}.inp".format(args.exe,args.single)
run=shlex.split(cmd)
with open("output.yml","wb") as out:
   subprocess.Popen(run,stdout=out)
# Banner
print("""\
========================================
Summary of the Test Suite results
========================================

Introduction
------------

A test is successful only if its result is "PASS". 
If one test is marked as Fail all test suite Failed


Test suite parameters
---------------------

The test comparison script has been run with the following parameters:

- Workdir       : %s
- Configuration : %s
- Reference dir : %s
- Test dir      : %s
- Single test   : %s
- Output file   : %s


Results
-------
"""%(os.getcwd(), args.config, args.refdir, args.testdir, args.single, args.output))

# Load global configuration
yaml_doc = YAML()
yaml_cfg = {"tests": []}
with open(args.config, "r") as cfg_file:
    yaml_cfg = yaml_doc.load(cfg_file)

# Perform test comparisons
siesta_tests = []
siesta_index = {}
siesta_xfail = []
for tcase in yaml_cfg["tests"]:
  if ("%s.yml"%(tcase["name"])==args.single): 
        # Init
    siesta_index[tcase["name"]] = len(siesta_tests)
    siesta_tests.append({"name": tcase["name"], "title": tcase["title"]})
    if ( ("fail_expect" in tcase) and (tcase["fail_expect"] == "yes") ):
        siesta_xfail.append(tcase["name"])

    # Check that test reference exists
    ref_path = os.path.join(args.refdir, "%s.yml" % tcase["name"])
    if ( not os.path.exists(ref_path) ):
        siesta_tests[-1]["result"] = "eref"
        siesta_tests[-1]["message"] = "missing reference file"
        continue

    # Load test reference
    yaml_doc = YAML()
    with open(ref_path, "r") as ref_file:
        yaml_ref = yaml_doc.load(ref_file)
    # Check that test output exists
    out_path = os.path.join(args.testdir,tcase["name"], args.yaml)
    if (args.single != "."):
        out_path = os.path.join(args.testdir,args.single)
    #out_path = os.path.join(args.testdir,tcase["name"]+".out")
    if ( not os.path.exists(out_path) ):
        siesta_tests[-1]["result"] = "eout"
        siesta_tests[-1]["message"] = "missing output file"
        continue

    # Load test output
    yaml_doc = YAML()
    with open(out_path, "r") as out_file:
        yaml_out = yaml_doc.load(out_file)

    # Compare energies
    siesta_tests[-1]["tolerances"] = {}
    ref_vars = yaml_ref["energies"]
    out_vars = yaml_out["energies"]
    tc_good = []
    tc_fail = []
    tc_skip = []
    print("{0:<s}: {1:s}-{2:s}\n".format("Detailed results for test",tcase["name"],tcase['title']))
    print("#{0:>10s} {1:>16s} {2:>16s} {3:>16s} {4:>16s} [{5:>6s}]".format("Observable","Expected","Actual","Difference","Tolerance","Status"))
    for (key, val) in ref_vars.items():
        if ( ("tolerances" in tcase.keys()) and (key in tcase["tolerances"]) ):
            tc_tol = tcase["tolerances"][key]
        elif ( key in yaml_cfg["tolerances"] ):
            tc_tol = yaml_cfg["tolerances"][key]
        else:
            tc_tol = args.abs_tol
        siesta_tests[-1]["tolerances"][key] = tc_tol
        isok=""
        if ( (key in out_vars) and (abs(out_vars[key] - val) < tc_tol) ):
            tc_good.append(key)
            isok="\033[92m Pass \033[0m"
        else:
            tc_fail.append(
                    {"name": key, "value": out_vars[key], "expected": val})
            isok="\033[91m Fail \033[0m"
        print("#{0:>10s} {1:16.8f} {2:16.8f} {3:16.8f} {4:16.8f} [{5:>6s}]".format(key,val,out_vars[key],out_vars[key]-val,tc_tol,isok))
    
            # Store comparison results
    tc_skip = sorted([item for item in out_vars.keys() if not item in ref_vars.keys()])

    # Store final test result
    if ( len(tc_fail) == 0 ):
        print("# Status: PASSED")
    else:
        print("# Status: FAILED")
    break

# Display statistics
ntests = len(tc_fail)+len(tc_skip)+len(tc_good)
ntpass = len(tc_good)
ntfail = len(tc_fail)
ntskip = len(tc_skip)
print("""\


Statistics
----------

Global results of the test suite:

- number of tests  : %4d
- successful tests : %4d (%5.1f%%)
- failed tests     : %4d (%5.1f%%)
- skiped tests     : %4d (%5.1f%%)

""" % (ntests, ntpass, ntpass*100.0/ntests, ntfail, ntfail*100.0/ntests,ntskip,ntskip*100/ntests))

# Display footer
print("""\
        .. note::

   This report is valid ReStructuredText, unless you ask for colorized output.
   You can render it in HTML, PDF and other formats using the Docutils package.
   Please consult `http://docutils.sourceforge.net/`_ for details.

Document generated on %s (UTC).
""" % strftime("%Y/%m/%d %H:%M:%S +0000",gmtime()))
