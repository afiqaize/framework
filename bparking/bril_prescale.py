#! /bin/env python
# usable only where brilws python3 env is valid
# inquires bril about trigger prescales ala 'brilcalc trg'
# https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html

import os
import sys
import json
import csv

from collections import OrderedDict
from argparse import ArgumentParser

def syscall(cmd, verbose = True):
    if verbose:
        print ("Executing: %s" % cmd)
    retval = os.system(cmd)
    if retval != 0:
        raise RuntimeError("Command failed!")

def rm_hlt_version(name):
    version_start = name.rfind("_v")
    if version_start == -1: 
        return name
    else:
        return name[:version_start+2]  

def parse_grl(grl):
        with open(grl) as ff:
                result = json.load(ff)

        runs = []
        for rr, ll in result.items():
                runs.append(rr)
        return sorted(runs)

def parse_bril(csvf, prescales):
        with open(csvf) as ff:
                pp = csv.DictReader(ff)

                for p in pp:
                        run = int(p["# run"])
                        ls = p["cmsls"]
                        if ls == "None":
                                continue
                        else:
                            ls = int(ls)

                        l1l = p["logic"]
                        hltnp = p["hltpath/prescval"].split('/')
                        hltn = rm_hlt_version(hltnp[0])
                        hltp = int(hltnp[1])

                        if l1l != "OR" and l1l != "ONE":
                                print("run", run, "ls", ls, "path", hltn, "l1 logic is", l1l, "instead of OR/ONE! unhandled case, skipping...")
                                continue

                        if hltp == 0:
                                continue

                        l1sp = p["l1bit/prescval"].split(' ')

                        l1s = []
                        l1p = []
                        for ss in l1sp:
                                sp = ss.split('/')

                                if int(sp[1]) != 0:
                                        l1s.append(sp[0])
                                        l1p.append(int(sp[1]))

                        if len(l1p) == 0:
                                continue

                        if run not in prescales:
                                prescales[run] = OrderedDict()

                        if ls not in prescales[run]:
                                prescales[run][ls] = OrderedDict()

                        if hltn in prescales[run][ls]:
                                print("run", run, "ls", ls, "path", hltn, "already in the prescale map! skipping...")
                                continue

                        prescales[run][ls][hltn] = OrderedDict()
                        prescales[run][ls][hltn]["hlt_prescale"] = hltp
                        prescales[run][ls][hltn]["seeds"] = l1s
                        prescales[run][ls][hltn]["seed_prescales"] = l1p

def inquire_bril(prescales, runs, triggers):
        nr = len(runs)
        for ir, run in enumerate(runs):
                if nr < 10 or ir % int(nr / 10) == 0:
                        print("processing", ir, "/", nr, "run...")

                for trigger in triggers:
                        syscall("brilcalc trg -r {run} --hltpath '{trigger}' --prescale -o test.csv".format(
                                run = run,
                                trigger = trigger
                        ), False)
                        parse_bril("test.csv", prescales)
                        syscall("rm test.csv", False)

if __name__ == '__main__':
        parser = ArgumentParser()
        parser.add_argument('--triggers', help = 'trigger expressions, comma separated', required = True)
        parser.add_argument('--runs', help = 'runs to consider, comma separated', default = '', required = False)
        parser.add_argument('--grl', help = 'good runs list e.g. golden json', default = '', required = False)
        parser.add_argument('--output', help = '', default = 'prescales.json', required = False)
        args = parser.parse_args()

        if (args.runs == '' and args.grl == '') or (args.runs != '' and args.grl != ''):
                raise RuntimeError("only one of --runs and --grl should be provided!")

        runs = sorted(args.runs.strip().split(',')) if args.runs != '' else parse_grl(args.grl)
        triggers = sorted(args.triggers.strip().split(','))

        prescales = OrderedDict()
        inquire_bril(prescales, runs, triggers)
        with open(args.output, "w") as jj: 
                json.dump(prescales, jj, indent = 1)
