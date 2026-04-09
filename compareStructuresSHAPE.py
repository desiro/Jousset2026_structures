#!/usr/bin/env python3
# script: compareStructuresSHAPE.py
# author: Daniel Desiro'
script_usage="""
usage
    compareStructuresSHAPE.py -pfx <out_prefix> -fs1 <fasta_genome_1> -fs2 <fasta_genome_2> -ds1 <shape_dir_1> -ds2 <shape_dir_2> -sf1 <structure_file_1> -sf2 <structure_file_2> [options]

version
    compareStructuresSHAPE.py 0.0.1

dependencies
    Python v3.9.7, ViennaRNA v2.5.0

################################################################

--prefix,-pfx
    output directory and prefix for result files

--genomeFasta1,-fs1
    first input genome FASTA file

--genomeFasta2,-fs2
    second input genome FASTA file

--dataSHAPE1,-sh1
    first directory with all SHAPE-MaP files; names have to match the fasta headers

--dataSHAPE2,-sh2
    second directory with all SHAPE-MaP files; names have to match the fasta headers

--structureFile1,-sf1
    first pickled file from getStructuresSHAPE

--structureFile2,-sf2
    second pickled file from getStructuresSHAPE

################################################################################

--columnSHAPE,-csh
    define column in SHAPE-MaP file (default: Norm_profile)

--reverse,-rev
    creates reverse of each sequence in the genome file (default: False)

--complement,-cmp
    creates complements of each sequence in the genome file (default: False)

################################################################################

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--pickleData,-pcl
    load from pickle file if available for faster setting changes 
    (default: False)

--junctionHeader,-jch
    define the junction output header and table data columns
    (default: aSeq,ai,aj,bSeq,bi,bj,dG,aMeanSHAPE,bMeanSHAPE,kmer,abDist,aDist,bDist,aSHAPE,bSHAPE,RNA,pattern,aLen,bLen)

--plotDistributions,-pdh
    plot junction distributions (default: False)

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
from matplotlib import pyplot as plt
from math import isnan

################################################################################
## main
################################################################################

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
    fasta_dict_1 = readFasta(opt, opt.fs1)
    fasta_dict_2 = readFasta(opt, opt.fs2)
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## read SHAPE-MaP files
    text_s = printLog(f"Status: Read SHAPE-MaP files ...", opt)
    shape_dict_1 = readSHAPE(opt, opt.sh1, fasta_dict_1)
    shape_dict_2 = readSHAPE(opt, opt.sh2, fasta_dict_2)
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## read structure files
    text_s = printLog(f"Status: Read structure SHAPE-MaP files ...", opt)
    (junction_list_1, junction_list_2), mutants = pickleParsing(f"raw", readStructureSHAPE, opt)
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## compare SHAPE structure
    text_s = printLog(f"Status: Compare structure SHAPE-MaP files ...", opt)
    unique_list_1, unique_list_2, overlap_list = pickleParsing(f"compare", compareStructureSHAPE, opt, *[junction_list_1, junction_list_2])
    time_s = getTime(opt, time_s, text_s)
    ########################################################################
    ## finalize junctions
    text_s = printLog(f"Status: Finalize and write SHAPE-MaP data ...", opt)
    ext1, ext2 = os.path.splitext(os.path.basename(os.path.abspath(opt.fs1)))[0], os.path.splitext(os.path.basename(os.path.abspath(opt.fs2)))[0]
    jcu_list_1 = (f"{ext1}_unique", unique_list_1, fasta_dict_1, shape_dict_1)
    jcu_list_2 = (f"{ext2}_unique", unique_list_2, fasta_dict_2, shape_dict_2)
    jco_list_1 = (f"{ext1}_overlap", overlap_list, fasta_dict_1, shape_dict_1)
    jco_list_2 = (f"{ext2}_overlap", overlap_list, fasta_dict_2, shape_dict_2)
    alist_1 = unique_list_1+overlap_list
    alist_2 = unique_list_2+overlap_list
    jca_list_1 = (f"{ext1}", alist_1, fasta_dict_1, shape_dict_1)
    jca_list_2 = (f"{ext2}", alist_2, fasta_dict_2, shape_dict_2)
    junction_dict = dict()
    for ext,junction_list,fasta_dict,shape_dict in [jcu_list_1, jcu_list_2, jco_list_1, jco_list_2, jca_list_1, jca_list_2]:
        junction_list = pickleParsing(f"finalized_{ext}", finalizeJunctionsSHAPE, opt, *[junction_list, fasta_dict, shape_dict, mutants])
        junction_dict[ext] = junction_list
        writeExtractions(opt, junction_list, ext)
        plotDotplotMulti(opt, [jc.dG for jc in junction_list], f"dG", f"distribution_{ext}_all")
        plotDotplotMulti(opt, [(jc.aMeanSHAPE+jc.aMeanSHAPE)/2 for jc in junction_list], f"mean SHAPE", f"distribution_{ext}_all")
    all_junctions_WT = junction_dict[f"{ext1}_unique"] + junction_dict[f"{ext1}_overlap"]
    plotDotplotMulti(opt, [jc.dG for jc in all_junctions_WT], f"dG", f"distribution_{ext1}")
    plotDotplotMulti(opt, [(jc.aMeanSHAPE+jc.aMeanSHAPE)/2 for jc in all_junctions_WT], f"mean SHAPE", f"distribution_{ext1}")
    all_junctions_MT = junction_dict[f"{ext2}_unique"] + junction_dict[f"{ext2}_overlap"]
    plotDotplotMulti(opt, [jc.dG for jc in all_junctions_MT], f"dG", f"distribution_{ext2}")
    plotDotplotMulti(opt, [(jc.aMeanSHAPE+jc.aMeanSHAPE)/2 for jc in all_junctions_MT], f"mean SHAPE", f"distribution_{ext2}")
    time_s = getTime(opt, time_s, text_s)
    ############################################################################
    return opt.pfx

################################################################################
## read fasta file

def readFasta(opt, fasta_file):
    ## read fasta file
    fasta_dict, RNA = dict(), ""
    with open(fasta_file, "r") as infa:
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

def readSHAPE(opt, shape_dir, fasta_dict):
    ## read SHAPE-MaP files
    shape_dict = dict()
    for name in fasta_dict.keys():
        shape_list = list()
        with open(os.path.join(shape_dir, f"{name}.txt"), "r") as insh:
            col = next(insh).strip().split().index(opt.csh)
            for line in insh:
                line = float(line.strip().split()[col])
                if isnan(line): line = 0.0
                shape_list.append(line)
        shape_dict[name] = shape_list
    return shape_dict

################################################################################
## read structure files

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

def readStructureSHAPE(opt):
    ## read structure files
    opt_pfx = opt.pfx
    data_list = list()
    mutants = list()
    for file in [opt.sf1, opt.sf2]:
        opt.pfx, fname = os.path.split(os.path.abspath(file))
        pickleFile, ext = os.path.splitext(fname)
        junction_list = loadData(pickleFile, opt)
        for jc in junction_list:
            if len(jc.aSeq.split("-")) > 1:
                jc.aSeq, mType = jc.aSeq.split("-")
                mutants.append((jc.aSeq, mType))
            else: jc.aSeq = jc.aSeq.split("-")[0]
            if len(jc.bSeq.split("-")) > 1:
                jc.bSeq, mType = jc.bSeq.split("-")
                mutants.append((jc.bSeq, mType))
            else: jc.bSeq = jc.bSeq.split("-")[0]
        data_list.append(junction_list)
    opt.pfx = opt_pfx
    return data_list, list(set(mutants))

################################################################################
## compare SHAPE structure

def compareStructureSHAPE(opt, junction_list_1, junction_list_2):
    ## compare structure files
    unique_list_1, unique_list_2, overlap_list = list(), list(), list()
    kmer_list = sorted(list(set([jc1.kmer for jc1 in junction_list_1]+[jc2.kmer for jc2 in junction_list_2])))
    for kmer in kmer_list:
        jc_list_1 = [jc1 for jc1 in junction_list_1 if jc1.kmer == kmer]
        jc_list_2 = [jc2 for jc2 in junction_list_2 if jc2.kmer == kmer]
        total1, total2 = len(jc_list_1), len(jc_list_2)
        for it,jc1 in enumerate(jc_list_1):
            if int(it*10000/total1) % 10 == 0 and int(it*100/total1) != 0:
                print(f"Status: Comparing {kmer:>2d} jc1 ... {it*100/total1:>4.1f} %             ", end="\r")
            if testCompare(jc1, jc_list_2): unique_list_1.append(jc1)
            else:                           overlap_list.append(jc1)
        for it,jc2 in enumerate(jc_list_2):
            if int(it*10000/total2) % 10 == 0 and int(it*100/total2) != 0:
                print(f"Status: Comparing {kmer:>2d} jc2 ... {it*100/total2:>4.1f} %             ", end="\r")
            if testCompare(jc2, jc_list_1): unique_list_2.append(jc2)
    return unique_list_1, unique_list_2, overlap_list

def testCompare(jcA, jc_list_B):
    ## co,pare junctions
    for jcB in jc_list_B:
        if   jcA.aSeq != jcB.aSeq: continue
        elif jcA.bSeq != jcB.bSeq: continue
        elif jcA.ai   != jcB.ai:   continue
        elif jcA.bi   != jcB.bi:   continue
        elif jcA.dG   != jcB.dG:   continue
        else: return False
    return True

################################################################################
## finalize junctions

def finalizeJunctionsSHAPE(opt, junction_list, fasta_dict, shape_dict, mutants):
    ## add shape and j positions and structure to junctions
    notWT = [k for k in fasta_dict.keys() if len(k.split("-")) > 1]
    for jc in junction_list:
        if notWT:
            for mutN,mutT in mutants:
                if jc.aSeq == mutN: jc.aSeq = f"{jc.aSeq}-{mutT}"
                if jc.bSeq == mutN: jc.bSeq = f"{jc.bSeq}-{mutT}"
        aSeq, bSeq = fasta_dict[jc.aSeq], fasta_dict[jc.bSeq]
        aShp, bShp = shape_dict[jc.aSeq], shape_dict[jc.bSeq]
        jc.aj, jc.bj = jc.ai+jc.kmer, jc.bi+jc.kmer
        jc.aSHAPE, jc.bSHAPE = aShp[jc.ai:jc.aj], bShp[jc.bi:jc.bj]
        jc.aMeanSHAPE, jc.bMeanSHAPE = mean(jc.aSHAPE), mean(jc.bSHAPE)
        jc.RNA = f"{aSeq[jc.ai:jc.aj]}&{bSeq[jc.bi:jc.bj]}"
        jc.pattern = "("*jc.kmer + "&" + ")"*jc.kmer
        jc.jType = "inter"
    return junction_list

def writeExtractions(opt, final_junctions, ext):
    ## write junctions
    out_name = os.path.splitext(os.path.basename(os.path.abspath(opt.pfx)))[0]
    out_write = os.path.join(opt.pfx, f"{out_name}_structures_inter{ext}.tsv")
    with open(out_write, "w") as outfile_inter:
        header = opt.jch.split(",")
        outfile_inter.write("\t".join(header)+"\n")
        for jc in sorted(final_junctions, key=attrgetter("dG","aSeq","bSeq","ai","bi","aj","bj"), reverse=False):
            line_list = list()
            for h in header:
                line, item = "", jc.__dict__[h]
                if   h in ["ai","bi"]:                        line = f"{item+1}"
                elif h in ["dG", "aMeanSHAPE", "bMeanSHAPE"]: line = f"{item:>.2f}"
                elif h in ["aSHAPE", "bSHAPE"]:               line = ";".join([f"{im:>.2f}" for im in item])
                else:                                         line = f"{item}"
                line_list.append(line)
            outfile_inter.write("\t".join(line_list)+"\n")

def plotDotplotMulti(opt, data_list, data, name):
    ## pyplot matrix plot
    outfile = os.path.join(opt.pfx, f"{name}_{data}")
    pdata = sorted(data_list)
    if data == "dG":
        yl = (-22,1)
    elif data == "mean SHAPE":
        yl = (-0.1,0.8)
    dv = std(pdata)
    mn = mean(pdata)
    lower = [mn-1*dv]*len(pdata)
    mnline = [mn]*len(pdata)
    upper = [mn+1*dv]*len(pdata)
    s = 1.0
    pz = plt.figure(figsize=(12*s,4*s))
    plt.plot(pdata, color="#d53e4f", linewidth=1, label=f"{data}")
    plt.plot(lower, color="#da9cf7", linestyle='--', linewidth=1, label="$\mu-\sigma$")
    plt.plot(mnline, color="#05aff7", linestyle='--', linewidth=1, label="$\mu$")
    plt.plot(upper, color="#da9cf7", linestyle=':', linewidth=1, label="$\mu+\sigma$")
    plt.ylim(yl)
    plt.xlabel("interactions")
    plt.ylabel(data)
    plt.legend()
    plt.subplots_adjust(bottom=0.2, top=1.2)
    pz.savefig(f"{outfile}.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{outfile}.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
    plt.close(pz)

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
