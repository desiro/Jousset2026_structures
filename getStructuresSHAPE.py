#!/usr/bin/env python3
# script: getStructuresSHAPE.py
# author: Daniel Desiro'
script_usage="""
usage
    getStructuresSHAPE.py -pfx <out_prefix> -fsa <fasta_genome> -dsh <shape_dir> [options]

version
    getStructuresSHAPE.py 0.0.1

dependencies
    Python v3.9.7, ViennaRNA v2.5.0

description
    This pipeline provides a framework for identifying short, structured, and thermodynamically
    stable RNA-RNA interactions using SHAPE-guided filtering and constrained folding. It further
    enables comparison of interaction landscapes across different datasets. Its main purpose was
    to detect stable continuous interactions for Jousset2026.

################################################################

--prefix,-pfx
    output directory and prefix for result files

--genomeFasta,-fsa
    input genome FASTA file

--dataSHAPE,-dsh
    directory with all SHAPE-MaP files; names have to match the fasta headers

################################################################################

--columnSHAPE,-csh
    define column in SHAPE-MaP file (default: Norm_profile)

--maxSHAPE,-xsh
    maximum allowed SHAPE-reactivity (default: 0.8)

--minRNA,-nra
    minimum allowed RNA interaction length (default: 4)

--maxRNA,-xra
    maximum allowed RNA interaction length (default: 20)

--interactionDistance,-itd
    maximum distance of two interacting regions from the UTR regions
    (default: 200)

--reverse,-rev
    creates reverse of each sequence in the genome file (default: False)

--complement,-cmp
    creates complements of each sequence in the genome file (default: False)

--structureEnergy,-sce
    minimum free energy threshold for a chimeric read structure in kcal/mol
    (default: -10.0)

################################################################################

--viennaRNAdangles,-vrd
    viennaRNA dangling ends setting for RNA folding (default: 2) (choices: 0,1,2,3)

--viennaRNAtemperature,-vrt
    viennaRNA temperature setting for RNA folding (default: 37.0)

--viennaRNAnoLP,-vrn
    viennaRNA no lonely pair setting for RNA folding (default: True)

################################################################################

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--pickleData,-pcl
    load from pickle file if available for faster setting changes 
    (default: False)

################################################################################

reference
    tba
"""

import sys
import os
import re
import time
import argparse as ap
import pickle
import traceback
from operator import attrgetter
from math import isnan
from itertools import product

try:
    from RNA import fold_compound, cvar, CONSTRAINT_DB, CONSTRAINT_DB_DEFAULT
except:
    print("Error: The mandatory library ViennaRNA 2.5 is missing, please read the README.md for more information!")
    exit()

################################################################################
## main

def main(opt):
    ############################################################################
    ## create output folder
    opt = makeDir(opt)
    time_s = getTime(opt)
    text_s = printLog(f"Status: Create output directory ...", opt)
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## read fasta file
    text_s = printLog(f"Status: Read fasta file ...", opt)
    fasta_dict = readFasta(opt)
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## read SHAPE-MaP files
    text_s = printLog(f"Status: Read SHAPE-MaP files ...", opt)
    shape_dict = readSHAPE(opt, fasta_dict)
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## get RNA k-mer
    text_s = printLog(f"Status: Get k-mers ...", opt)
    kmer_dict = pickleParsing(f"kmer", getKmer, opt, *[fasta_dict, shape_dict])
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## find rev-comp k-mers
    text_s = printLog(f"Status: Find reverse k-mers ...", opt)
    junction_list = pickleParsing(f"reverse", findReverseKmers, opt, *[kmer_dict, fasta_dict])
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## recalculate delta G
    text_s = printLog(f"Status: Recalculate dG ...", opt)
    junction_list = pickleParsing(f"recalculate", recalculateDeltaG, opt, *[junction_list, fasta_dict])
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## remove subsets
    text_s = printLog(f"Status: Remove subsets ...", opt)
    junction_list = pickleParsing(f"subsets", removeSubsets, opt, *[junction_list, fasta_dict])
    time_s = getTime(opt, time_s, text_s)
    ############################################################################
    return opt.pfx

################################################################################
## read fasta file

def readFasta(opt):
    ## read fasta file
    fasta_dict, RNA = dict(), ""
    with open(opt.fsa, "r") as infa:
        for line in infa:
            line = line.strip()
            if re.match(r">", line):
                if RNA: fasta_dict[name] = revComp(RNA, opt)
                RNA, name = "", line[1:].split()[0]
            else:
                RNA += line
        fasta_dict[name] = revComp(RNA, opt)
    return fasta_dict

def revComp(RNA, opt):
    ## complement dictionary, or transform DNA to RNA
    RNA = RNA.upper()
    D2Rc = {"A":"U","T":"A","U":"A","C":"G","G":"C","R":"Y","Y":"R","M":"K",\
            "K":"M","S":"W","W":"S","B":"V","V":"B","D":"H","H":"D","N":"N"}
    if opt.cmp: RNA = "".join(D2Rc[i] for i in RNA)
    else:              RNA = RNA.replace("T","U")
    if opt.rev: RNA = RNA[::-1]
    return RNA

################################################################################
## read SHAPE-MaP files

def readSHAPE(opt, fasta_dict):
    ## read SHAPE-MaP files
    shape_dict = dict()
    for name in fasta_dict.keys():
        shape_list = list()
        with open(os.path.join(opt.dsh, f"{name}.txt"), "r") as insh:
            col = next(insh).strip().split().index(opt.csh)
            for line in insh:
                line = float(line.strip().split()[col])
                if isnan(line): line = 0.0
                shape_list.append(line)
        shape_dict[name] = shape_list
    return shape_dict

################################################################################
## get RNA k-mer

def getKmer(opt, fasta_dict, shape_dict):
    ## get RNA k-mers
    kmer_dict = dict()
    for name,RNA in fasta_dict.items():
        xmer_dict = dict()
        for xmer in range(opt.nra, opt.xra+1):
            km_list = [testSHAPE(opt, RNA[i:i+xmer], i, xmer, shape_dict[name][i:i+xmer]) for i in range(0,len(RNA)-xmer+1)]
            km_list = [km for km in km_list if km is not None]
            km_dict = sortKmer(opt, km_list)
            xmer_dict[xmer] = km_dict
        kmer_dict[name] = xmer_dict
    return kmer_dict

def testSHAPE(opt, RNA, iv, xmer, SHAPE):
    ## test SHAPE
    if not len([0 for i in SHAPE if i > opt.xsh]): return((RNA, iv))

def sortKmer(opt, km_list):
    # sort kmers
    km_dict = dict()
    for RNA,i in km_list:
        km_dict[RNA] = km_dict.get(RNA, list()) + [i]
    return km_dict

################################################################################
## find rev-comp k-mers

def findReverseKmers(opt, kmer_dict, fasta_dict):
    ## find reverse complemetn k-mer
    order_dict = {k:i for i,k in enumerate(fasta_dict.keys())}
    multi_list = [(aName, bName, aRNA, bRNA, ailist, bilist, xmer) for aName,aXmer_dict in kmer_dict.items() for bName,bXmer_dict in kmer_dict.items() for xmer,aKm_dict in aXmer_dict.items() for aRNA,ailist in aKm_dict.items() for bRNA,bilist in bXmer_dict[xmer].items() if order_dict[aName] < order_dict[bName]]
    data_list = findReverseKmersMulti(multi_list, opt, 0, len(multi_list), list())
    return data_list

def findReverseKmersMulti(multi_list, opt, run, split, res_list):
    ## find reverse complemetn k-mer
    for it,mltlst in enumerate(multi_list):
        if int(it*10000/split) % 10 == 0 and int(it*100/split) != 0:
            print(f"Status: Folding {run:>2d} ... {it*100/split:>4.1f} %             ", end="\r")
        aName, bName, aRNA, bRNA, ailist, bilist, xmer = mltlst
        n = "n"*opt.vrd
        RNA = f"{n}{aRNA}{n}&{n}{bRNA}{n}"
        dG, pattern = doCofold(RNA, f"..{'<'*len(aRNA)}....{'>'*len(bRNA)}..", opt)
        if len(aRNA) == pattern.count("("):
            RNA = f"{aRNA}&{bRNA}"
            pattern = pattern[:len(aRNA.split("&")[0])+(2*opt.vrd)][opt.vrd:-opt.vrd]+"&"+pattern[len(bRNA.split("&")[0])+(2*opt.vrd):][opt.vrd:-opt.vrd]
            km = {"aSeq":aName, "ai":ailist, "bSeq":bName, "bi":bilist, "kmer":xmer, "dG":round(dG,2)}
            kmj = junction(**km)
            res_list.append(kmj)
    return res_list

def doCofold(RNA, constraint, opt):
    ## do Cofold
    cvar.dangles = opt.vrd
    cvar.noLP = int(opt.vrn)
    cvar.temperature = opt.vrt
    fc = fold_compound(RNA)
    fc.constraints_add(constraint, CONSTRAINT_DB | CONSTRAINT_DB_DEFAULT)
    pattern, mfe = fc.mfe()
    return mfe, pattern

################################################################################
## recalculate delta G

def recalculateDeltaG(opt, junction_list, fasta_dict):
    ## removes prefix and suffix
    junction_list = sorted(junction_list, key=attrgetter("kmer"))
    total = len(junction_list)
    new_junctions = list()
    flens = {k:len(v) for k,v in fasta_dict.items()}
    for it,jc in enumerate(junction_list):
        if int(it*10000/total) % 10 == 0 and int(it*100/total) != 0:
            print(f"Status: Folding ... {it*100/total:>4.1f} %             ", end="\r")
        aRNA, bRNA = fasta_dict[jc.aSeq][jc.ai[0]:jc.ai[0]+jc.kmer], fasta_dict[jc.bSeq][jc.bi[0]:jc.bi[0]+jc.kmer]
        n = "n"*opt.vrd
        RNA = f"{n}{aRNA}{n}&{n}{bRNA}{n}"
        dG, pattern = doPartition(RNA, f"..{'<'*len(aRNA)}....{'>'*len(bRNA)}..", opt)
        pattern = pattern[:len(aRNA.split("&")[0])+(2*opt.vrd)][opt.vrd:-opt.vrd]+"&"+pattern[len(bRNA.split("&")[0])+(2*opt.vrd):][opt.vrd:-opt.vrd]
        jc.dG = dG
        if jc.dG <= opt.sce:
            for ai,bi in product(jc.ai,jc.bi):
                distance, aDist, bDist, aLen, bLen = calculateDistance(jc, ai, bi, flens)
                if distance <= opt.itd:
                    km = {"aSeq":jc.aSeq, "ai":ai, "bSeq":jc.bSeq, "bi":bi, "kmer":jc.kmer, "dG":jc.dG, "abDist":distance, "aDist":aDist, "bDist":bDist, "aLen":aLen, "bLen":bLen}
                    kmj = junction(**km)
                    new_junctions.append(kmj)
    return new_junctions

class junction(object):
    def __init__(self, **data):
        self.__dict__.update((k,self.transf(v)) for k,v in data.items())
    def transf(self, s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    def plot(self, sep):
        ldat = sep.join([f"{var}" for key,var in vars(self).items()])
        return ldat
    def __eq__(self, other): 
        if not isinstance(other, junction):
            return NotImplemented
        return not [0 for sk,si in self.__dict__.items() if si != other.__dict__[sk]]

def doPartition(RNA, constraint, opt):
    ## do Cofold
    cvar.dangles = opt.vrd
    cvar.noLP = int(opt.vrn)
    cvar.temperature = opt.vrt
    cvar.pf_scale = -10.0
    fc = fold_compound(RNA)
    fc.constraints_add(constraint, CONSTRAINT_DB | CONSTRAINT_DB_DEFAULT)
    _, tmfe = fc.mfe()
    fc.exp_params_rescale(tmfe)
    fc.pf()
    bpp = fc.bpp()
    pattern, _, _, _, ensembleFreeEnergy = fc.pf_dimer()
    return ensembleFreeEnergy, pattern

def calculateDistance(jc, ai, bi, flens):
    # calculates distance of two interacting regions on the segment
    n2a, n2b = flens[jc.aSeq]/2, flens[jc.bSeq]/2
    adist = n2a - abs(n2a - (ai + jc.kmer/2))
    bdist = n2b - abs(n2b - (bi + jc.kmer/2))
    abdist = abs(adist-bdist)
    return abdist, adist, bdist, flens[jc.aSeq], flens[jc.bSeq]

################################################################################
## remove subsets

def removeSubsets(opt, junction_list, fasta_dict):
    ## removes prefix and suffix
    junction_list = sorted(junction_list, key=attrgetter("kmer"))
    junction_list_rev = sorted(junction_list, key=attrgetter("kmer"), reverse=True)
    total, jc_count, new_junctions = len(junction_list), 0, list()
    for it,jc in enumerate(junction_list):
        print(f"Status: junction {it}/{total} at kmer {jc.kmer} with extended {len(new_junctions)}/{jc_count}                 ", end="\r")
        jc_count += 1
        if testSuffixPrefix(jc, junction_list_rev):
            new_junctions.append(jc)
    print(f"Status: junctions {len(new_junctions)}/{jc_count}                                                  ")
    return new_junctions

def testSuffixPrefix(jc, junction_list_rev):
    ## test if jc is a suffix or prefix
    for jct in junction_list_rev:
        if jc.kmer >= jct.kmer: break
        if jc.aSeq != jct.aSeq or jc.bSeq != jct.bSeq: continue
        if jct.ai <= jc.ai and jc.ai+jc.kmer <= jct.ai+jct.kmer and jct.bi <= jc.bi and jc.bi+jc.kmer <= jct.bi+jct.kmer:
            return False
    return True

################################################################################
## basic functions

def makeDir(opt):
    ## create directory
    dir_name, dir_base = opt.pfx, opt.pfx
    if not opt.ovr:
        i = 1
        while os.path.isdir(dir_name):
            dir_name = f"{dir_base}_{i}"
            i += 1
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    opt.pfx = dir_name
    return opt

def getTime(opt, time_s=0, name=""):
    if time_s and name:
        time_e = time.time()-time_s
        time_e = time.strftime("%H:%M:%S", time.gmtime(time_e))
        time_c = time.strftime('%x %X')
        printLog(f"Status: {name} finished at {time_c} in {time_e}", opt)
    return time.time()

def printLog(status, opt):
    print(status)
    if opt.log:
        f_name = os.path.splitext(os.path.basename(os.path.abspath(opt.pfx)))[0]
        with open(os.path.join(opt.pfx, f"{f_name}.log"), "a") as logout:
            logout.write(f"{status}\n")
    return status.replace("Status: ","").replace(" ...","")

def pickleParsing(pickleFile, parseFunction, opt, *argv):
    ## parses pickle files and function calling
    pickleFile = f"{os.path.basename(os.path.abspath(opt.pfx))}_{pickleFile}"
    if opt.pcl and os.path.isfile(os.path.join(opt.pfx, f"{pickleFile}.pcl")):
        pickle_data = loadData(pickleFile, opt)
    else:
        pickle_data = parseFunction(opt, *argv)
        saveData(pickle_data, pickleFile, opt)
    return pickle_data

def loadData(fname, opt):
    ## load data with pickle
    pcl_file = os.path.join(opt.pfx, f"{fname}.pcl")
    with open(pcl_file, "r+b") as pcl_in:
        pcl_data = pickle.load(pcl_in)
    return pcl_data

def saveData(pcl_data, fname, opt):
    ## save data with pickle
    pcl_file = os.path.join(opt.pfx, f"{fname}.pcl")
    with open(pcl_file, "w+b") as pcl_out:
        pickle.dump(pcl_data, pcl_out , protocol=4)

################################################################################
## parser

class options(object):
    def __init__(self, **data):
        self.__dict__.update((k,v) for k,v in data.items())
    def plot(self, sep):
        ldat = sep.join([f"{var}" for key,var in vars(self).items()])
        return ldat

if __name__ == "__main__":

    ############################################################################
    ## get time and save call
    sscript = sys.argv[0]
    start_time = time.time()
    current_time = time.strftime('%x %X')
    scall = " ".join(sys.argv[1:])
    with open(f"{sscript}.log", "a") as calllog:
        calllog.write(f"Start : {current_time}\n")
        calllog.write(f"Script: {sscript}\n")
        calllog.write(f"Call  : {scall}\n")
    print(f"Call: {scall}")
    print(f"Status: Started at {current_time}")
    ############################################################################
    ## transform string into int, float, bool if possible
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    ############################################################################
    ## save documentation
    rx_text = re.compile(r"\n^(.+?)\n((?:.+\n)+)",re.MULTILINE)
    rx_oneline = re.compile(r"\n+")
    rx_options = re.compile(r"\((.+?)\:(.+?)\)")
    help_dict, type_dict, text_dict, mand_list = {}, {}, {}, []
    for match in rx_text.finditer(script_usage):
        argument = match.groups()[0].strip()
        text = " ".join(rx_oneline.sub("",match.groups()[1].strip()).split())
        argopts = {"action":"store", "help":None, "default":None, "choices":None}
        for option in rx_options.finditer(text):
            key = option.group(1).strip()
            var = option.group(2).strip()
            if var == "False": argopts["action"] = "store_true"
            if var == "True": argopts["action"] = "store_false"
            if key == "choices": var = [vs.strip() for vs in var.split(",")]
            if key == "default": var = trans(var)
            argopts[key] = var
        if argopts["default"]: add_default = f" (default: {str(argopts['default'])})"
        else: add_default = ""
        argopts["help"] = rx_options.sub("",text).strip()+add_default
        argnames = argument.split(",")
        if len(argnames) > 1:
            if argopts["default"] == None:
                mand_list.append(f"{argnames[1][1:]}")
            type_dict[f"{argnames[1][1:]}"] = argopts["default"]
            argopts["argshort"] = argnames[1]
            help_dict[argnames[0]] = argopts
        else:
            text_dict[argnames[0]] = argopts["help"]
    ############################################################################
    ## get arguments
    if text_dict["dependencies"]:
        desc = f"{text_dict['description']} (dependencies: {text_dict['dependencies']})"
    else:
        desc = text_dict['description']
    p = ap.ArgumentParser(prog=sscript, prefix_chars="-", usage=text_dict["usage"],
                          description=desc, epilog=text_dict["reference"])
    p.add_argument("-v", "--version", action="version", version=text_dict["version"])
    for argname,argopts in help_dict.items():
        argshort = argopts["argshort"]
        if argopts["choices"]:
            p.add_argument(argshort, argname,            dest=f"{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"],   choices=argopts["choices"])
        else:
            p.add_argument(argopts["argshort"], argname, dest=f"{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"])
    p._optionals.title = "arguments"
    opt = vars(p.parse_args())
    ############################################################################
    ## validate arguments
    if None in [opt[mand] for mand in mand_list]:
        print("Error: Mandatory arguments missing!")
        print(f"Usage: {text_dict['usage']} use -h or --help for more information.")
        sys.exit()
    for key,var in opt.items():
        if key not in mand_list:
            arg_req, arg_in = type_dict[key], trans(var)
            if type(arg_req) == type(arg_in):
                opt[key] = arg_in
            else:
                print(f"Error: Argument {key} is not of type {type(arg_req)}!")
                sys.exit()
    ############################################################################
    ## add log create options class
    opt["log"] = True
    copt = options(**opt)
    ############################################################################
    ## call main function
    try:
        #saved = main(opt)
        saved = main(copt)
    except KeyboardInterrupt:
        print("Error: Interrupted by user!")
        sys.exit()
    except SystemExit:
        print("Error: System exit!")
        sys.exit()
    except Exception:
        print("Error: Script exception!")
        traceback.print_exc(file=sys.stderr)
        sys.exit()
    ############################################################################
    ## finish
    started_time = current_time
    elapsed_time = time.time()-start_time
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    current_time = time.strftime('%x %X')
    if saved:
        with open(f"{sscript}.log", "a") as calllog,\
             open(os.path.join(saved,f"call.log"), "a") as dirlog:
            calllog.write(f"Save  : {os.path.abspath(saved)}\n")
            calllog.write(f"Finish: {current_time} in {elapsed_time}\n")
            ## dirlog
            dirlog.write(f"Start : {started_time}\n")
            dirlog.write(f"Script: {sscript}\n")
            dirlog.write(f"Call  : {scall}\n")
            dirlog.write(f"Save  : {os.path.abspath(saved)}\n")
            dirlog.write(f"Finish: {current_time} in {elapsed_time}\n")
    print(f"Status: Saved at {saved}")
    print(f"Status: Finished at {current_time} in {elapsed_time}")
    sys.exit(0)
