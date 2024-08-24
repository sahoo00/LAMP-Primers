import io
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import json
import re
from dna_features_viewer import GraphicFeature, GraphicRecord
import pandas as pd
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import re
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils import MeltingTemp as mt

def getGC(seq):
    try:
        from Bio.SeqUtils import GC
        return GC(seq)
    except ImportError:
        from Bio.SeqUtils import gc_fraction
        return gc_fraction(seq, "weighted")*100

def splitPrimers(target, primer):
    i1 = target.find(primer[-15:])
    i2 = len(primer) - 15
    while (i1 > 0 and i2 > 0 and target[i1] == primer[i2]):
        i1 -= 1
        i2 -= 1
    return [primer[:(i2+1)], primer[(i2+1):]]

def getPos(seq, primer):
    i1 = seq.find(primer)
    if i1 < 0:
        i1 = seq.find(Seq(str(primer)).reverse_complement())
    return i1

def binarySearch(array, target):
    left, right = 0, len(array)-1
    while left <= right:
        middle = (left+right) // 2
        potentialMatch = array[middle]
        if target == potentialMatch[1]:
            return middle
        elif target < potentialMatch[1]:
            right = middle - 1
        else:
            left = middle + 1
    return left
def getSeq(seqrecord):
    if seqrecord is None:
        return None
    else:
        return seqrecord.seq

def getTm(seq):
    res = mt.Tm_NN(seq, dnac1=30.14, dnac2=0,
                     Na=50, Mg=4.0, Tris=0, dNTPs=0,
                     nn_table=mt.DNA_NN3, saltcorr=1)
    return res
def CountNeighbors(theSeekSeq, theSequence):
    returnValue = 0;
    i = 0;
    while (i>=0 and i<len(theSequence)):
        try:
            i=theSequence.index(theSeekSeq,i);
            if (i>=0):
                returnValue+= 1;
                i+= 1;
        except Exception as inst:
            break
    return returnValue
def deltaG(seq):
    aaCount = CountNeighbors("AA",seq)+CountNeighbors("TT",seq)+CountNeighbors("UU",seq);
    atCount = CountNeighbors("AT",seq)+CountNeighbors("AU",seq);
    taCount = CountNeighbors("TA",seq)+CountNeighbors("UA",seq);
    caCount = CountNeighbors("CA",seq)+CountNeighbors("TG",seq)+CountNeighbors("UG",seq);
    gtCount = CountNeighbors("GT",seq)+CountNeighbors("AC",seq)+CountNeighbors("GU",seq);
    ctCount = CountNeighbors("CT",seq)+CountNeighbors("AG",seq)+CountNeighbors("CU",seq);
    gaCount = CountNeighbors("GA",seq)+CountNeighbors("TC",seq)+CountNeighbors("UC",seq);
    cgCount = CountNeighbors("CG",seq);
    gcCount = CountNeighbors("GC",seq);
    ggCount = CountNeighbors("GG",seq)+CountNeighbors("CC",seq);
    counts = np.array([aaCount, atCount, taCount, caCount, gtCount, ctCount,
             gaCount, cgCount, gcCount, ggCount])
    table = np.array([-1.00, -0.88, -0.58, -1.45, -1.44, -1.28,
                      -1.30, -2.17, -2.24, -1.84])
    initGC = 0.98; initAT = 1.03; sym = 0;
    #table = np.array([-1.02, -0.73, -0.60, -1.38, -1.43,
    #                  -1.16, -1.46, -2.09, -2.28, -1.77])
    #initGC = 1.82; initAT = 2.8; sym = 0.4;
    val = np.sum(counts*table)
    ahash = {'A', 'T', 'U'}
    chash = {'G', 'C'}
    if seq[-1] in ahash:
        val += initAT
    else:
        val += initGC
    if seq[0] in ahash:
        val += initAT
    else:
        val += initGC
    if str(Seq(str(seq)).reverse_complement()) == seq:
        val += sym
    return val
class LAMP:
    reactions = {}
    reactions["Na"] = 50.0
    reactions["Mg"] = 4.0
    lengths = {}
    lengths["F1c/B1c"] = [20, 22]
    lengths["F2/B2"] = [18, 20]
    lengths["F3/B3"] = [18, 20]
    lengths["LF/LB"] = [15, 25]
    Tm = {}
    Tm["F1c/B1c"] = [64, 66]
    Tm["F2/B2"] = [59, 61]
    Tm["F3/B3"] = [59, 61]
    Tm["LF/LB"] = [64, 66]
    GC = {}
    GC["F1c/B1c"] = [40, 65]
    GC["F2/B2"] = [40, 65]
    GC["F3/B3"] = [40, 65]
    GC["LF/LB"] = [40, 65]
    stability = {}
    stability["5p"] = -4
    stability["3p"] = -4
    stability["3p loop"] = -2.0
    stability["dimer"] = -2.5
    stability["dimer loop"] = -3.5
    distances = {}
    distances["F2-B2"] = [120, 160]
    distances["loop F1c-F2"] = [40, 60]
    distances["F2-F3"] = [0, 60]
    distances["F1c-B1c"] = [0, 100]
    limits = {}
    limits["F1c/B1c"] = 3
    limits["F2/B2"] = 10
    limits["F3/B3"] = 3
    limits["LF/LB"] = 10
    sets = 1000

    def initSet(self, pset):
        f3, f2, f1, b1, b2, b3 = pset
        self.pset = pset
        self.F3 = SeqRecord(f3[0], id="F3", name="F3", description="F3")
        self.F2 = SeqRecord(f2[0], id="F2", name="F2", description="F2")
        self.F1c = SeqRecord(f1[0].reverse_complement(), id="F1c", name="F1c", description="F1c")
        self.B2 = SeqRecord(b2[0].reverse_complement(), id="B2", name="B2", description="B2")
        self.B1c = SeqRecord(b1[0], id="B1c", name="B1c", description="B1c")
        self.B3 = SeqRecord(b3[0].reverse_complement(), id="B3", name="B3", description="B3")
        self.LB = None
        self.LF = None
        self.FIP = SeqRecord(self.getFIP(), id="FIP", name="FIP", description="FIP")
        self.BIP = SeqRecord(self.getBIP(), id="BIP", name="BIP", description="BIP")
        self.id = 0
        self.name = "LAMP"
        self.id1 = ""
        return
    def initText(self, text):
        handle = io.StringIO(text)
        primers = SeqIO.parse(handle, 'fasta')
        ahash = {s.id:s for s in primers}
        self.F3 = ahash['F3']
        self.B3 = ahash['B3']
        self.F2 = ahash['F2']
        self.F1c = ahash['F1c']
        self.B2 = ahash['B2']
        self.B1c = ahash['B1c']
        self.LB = None
        self.LF = None
        if 'LB' in ahash:
            self.LB = ahash['LB']
        if 'LF' in ahash:
            self.LF = ahash['LF']
        self.FIP = SeqRecord(self.getFIP(), id="FIP", name="FIP", description="FIP")
        self.BIP = SeqRecord(self.getBIP(), id="BIP", name="BIP", description="BIP")
        self.id = 0
        self.name = "LAMP"
        self.id1 = ""
        return
    def __init__(self, text):
        if type(text) == str:
            self.initText(text)
        else:
            self.initSet(text)
        return
    def addLF(self, lf):
        self.lfdata = lf
        self.LF = SeqRecord(lf[0].reverse_complement(), id="LF", name="LF", description="LF")
        return
    def addLB(self, lb):
        self.lbdata = lb
        self.LB = SeqRecord(lb[0], id="LB", name="LB", description="LB")
        return
    def getFIP(self):
        return self.F1c.seq + self.F2.seq
    def getBIP(self):
        return self.B1c.seq + self.B2.seq
    @staticmethod
    def getPTPRC_1():
        text = '''
>F3
ACAAAAACTTCCCCAGAAGA
>B3
CCTAGCTTTGCGTAGAGC
>FIP
GAGATCCATCCCTGCAGTGACTGAAGGGAACAAGCATCA
>BIP
CAAACGGGAATATTTTGTGCTTTGTTACCACTTGAAAAATATCCACTAC
>LB
ATCTCTTAGAAAGTGCGGAAACAGA
>F2
CTGAAGGGAACAAGCATCA
>F1c
GAGATCCATCCCTGCAGTGA
>B2
TACCACTTGAAAAATATCCACTAC
>B1c
CAAACGGGAATATTTTGTGCTTTGT
'''
        return LAMP(text)
    @staticmethod
    def getKCNJ15_1():
        text = '''
>F3
GCCAGGAGTCTGCACTCT
>B3
AGCAGTGTGCTTCACCAG
>LF
CCTGGCTACCATGGGATTCTG
>LB
GAGTGACCAGTGTTTCCAGAGC
>F2
TCAGTCTTTGCAGGCAGTAG
>F1c
ACGTCCTCGCTCCCCTTCAC
>B2
CCGATGTGAATGGCATCCAT
>B1c
AAGAAGACACCTGACCTGCGG
'''
        return LAMP(text)
    @staticmethod
    def getPTPRCTarget():
        Entrez.email = "dsahoo@ucsd.edu"  # Always tell NCBI who you are
        term = 'PTPRC[gene] AND human[orgn] AND mrna[filter]'
        handle = Entrez.esearch(db="nucleotide", term=term, rettype="uilist", retmode="json")
        obs = json.load(handle)
        id1 = obs["esearchresult"]["idlist"][2]
        handle = Entrez.efetch(db="nucleotide", id=id1, rettype="fasta", retmode="text")
        from Bio import SeqIO
        seq = list(SeqIO.parse(handle, 'fasta'))[0]
        return seq[3349:5000].seq
    def getKCNJ15Target():
        Entrez.email = "dsahoo@ucsd.edu"  # Always tell NCBI who you are
        #term = 'KCNJ15[gene] AND human[orgn] AND mrna[filter]'
        #handle = Entrez.esearch(db="nucleotide", term=term, rettype="uilist", retmode="json")
        #obs = json.load(handle)
        #id1 = obs["esearchresult"]["idlist"][7]
        id1 = "1519244435"
        handle = Entrez.efetch(db="nucleotide", id=id1, rettype="fasta", retmode="text")
        from Bio import SeqIO
        seq = list(SeqIO.parse(handle, 'fasta'))[0]
        return seq[:713].seq
    def getTYROBPTarget():
        Entrez.email = "dsahoo@ucsd.edu"  # Always tell NCBI who you are
        id1 = "NM_003332.4" #TYROBP
        handle = Entrez.efetch(db="nucleotide", id=id1, rettype="fasta", retmode="text")
        from Bio import SeqIO
        seq = list(SeqIO.parse(handle, 'fasta'))[0]
        return seq.seq
    def getC1QATarget():
        Entrez.email = "dsahoo@ucsd.edu"  # Always tell NCBI who you are
        id1 = "NM_015991.4" #C1QA
        handle = Entrez.efetch(db="nucleotide", id=id1, rettype="fasta", retmode="text")
        from Bio import SeqIO
        seq = list(SeqIO.parse(handle, 'fasta'))[0]
        return seq.seq
    def getC1qaTarget():
        Entrez.email = "dsahoo@ucsd.edu"  # Always tell NCBI who you are
        id1 = "NM_007572.2" #C1qa
        handle = Entrez.efetch(db="nucleotide", id=id1, rettype="fasta", retmode="text")
        from Bio import SeqIO
        seq = list(SeqIO.parse(handle, 'fasta'))[0]
        return seq.seq
    def amplifyTarget(self, target):
        targets = [target]
        return
    def __repr__(self):
        lkeys = ['F3', 'B3', 'F2', 'B2', 'F1c', 'B1c', 'LF', 'LB']
        id1 = self.id1
        name = self.name
        res = f"LAMP {id1} {name}\n"
        for k in lkeys:
            if self.__dict__[k] is not None:
                res += f"{k}={self.__dict__[k].seq}\n"
        return res
    def printParameters(self, seq):
        lkeys = ['F3', 'B3', 'F2', 'B2', 'F1c', 'B1c', 'LF', 'LB']
        chash = {'F3':0, 'B3':1, 'F2':0, 'B2':1, 'F1c':1, 'B1c':0, 'LF':1, 'LB':0}
        res = ''
        for k in lkeys:
            i1 = -1
            if self.__dict__[k] is not None:
                if chash[k] == 0:
                    i1 = seq.find(self.__dict__[k].seq)
                else:
                    i1 = seq.find(self.__dict__[k].seq.reverse_complement())
            if i1 < 0:
                continue
            res += f"{k} 5pos={i1} 3pos={i1+len(self.__dict__[k])-1} len={len(self.__dict__[k])}"
            res += f" Tm={getTm(self.__dict__[k].seq):0.2f} GC={getGC(self.__dict__[k].seq):.0f}"
            res += f" 5dg={deltaG(self.__dict__[k].seq[:6]):.2f}"
            res += f" 3dg={deltaG(self.__dict__[k].seq[-6:]):.2f}"
            res += "\n"
        print(res)
        return
    @staticmethod
    def printPrimerParameters(primer, seq=None):
        i1 = -1
        if seq is not None:
            i1 = seq.find(primer)
            if i1 < 0:
                i1 = seq.find(primer.reverse_complement())
        res = ''
        res += f"{primer} 5pos={i1} 3pos={i1+len(primer)-1} len={len(primer)}"
        res += f" Tm={getTm(primer):0.2f} GC={getGC(primer):.0f}"
        res += f" 5dg={deltaG(primer[:6]):.2f}"
        res += f" 3dg={deltaG(primer[-6:]):.2f}"
        res += "\n"
        print(res)
        return
    @staticmethod
    def getSeqParamaters(seq):
        seq = Seq(str(seq))
        seq.Tm = getTm(seq)
        seq.GC = getGC(seq)
        seq.dg5 = deltaG(seq[:6])
        seq.dg3 = deltaG(seq[-6:])
        return seq
    @staticmethod
    def isStable(seq):
        res = 1
        if seq.dg5 > __class__.stability["5p"]:
            res = 0
        if seq.dg3 > __class__.stability["3p"]:
            res = 0
        return res
    @staticmethod
    def isValid(seq, k):
        res = 1
        if len(seq) < __class__.lengths[k][0]:
            res = 0
        if len(seq) > __class__.lengths[k][1]:
            res = 0
        if seq.Tm < __class__.Tm[k][0]:
            res = 0
        if seq.Tm > __class__.Tm[k][1]:
            res = 0
        if seq.GC < __class__.GC[k][0]:
            res = 0
        if seq.GC > __class__.GC[k][1]:
            res = 0
        return res
    @staticmethod
    def generatePrimers(seq):
        minlen = __class__.distances["F2-B2"][0] + 2 * __class__.lengths["F3/B3"][0]
        if len(seq) < minlen:
            return None
        m1 = np.min(np.array([__class__.lengths[k] for k in __class__.lengths]))
        m2 = np.max(np.array([__class__.lengths[k] for k in __class__.lengths]))
        n = len(seq)
        primers = []
        f3b3s = []
        f2b2s = []
        f1b1s = []
        lflbs = []
        for pos in range(n - m1):
            for j in range(m1, m2+1):
                primer = __class__.getSeqParamaters(seq[pos:(pos+j)])
                if __class__.isStable(primer):
                    primers += [[primer, pos, j]]
                    if __class__.isValid(primer, "F3/B3"):
                        f3b3s += [[primer, pos, j]]
                    if __class__.isValid(primer, "F2/B2"):
                        f2b2s += [[primer, pos, j]]
                    if __class__.isValid(primer, "F1c/B1c"):
                        f1b1s += [[primer, pos, j]]
                    if __class__.isValid(primer, "LF/LB"):
                        lflbs += [[primer, pos, j]]
        res = []
        phash = {'F3':{},'F2':{},'LF':{},'F1':{},'B1':{},'LB':{},'B2':{},'B3':{}}
        f3b3sel = []
        for f3 in f3b3s:
            if f3[1] > (n - minlen):
                continue
            f3b3sel += [f3]
        f3b3sel = sorted(f3b3sel, key=lambda x: max(x[0].dg5, x[0].dg3))
        for f3 in f3b3sel:
            if f3[1] in phash['F3']:
                continue
            phash['F3'][f3[1]] = f3
            def getPrimerSet(f3):
                f2_s = binarySearch(f2b2s, f3[1] + f3[2] + __class__.distances["F2-F3"][0])
                f2_e = binarySearch(f2b2s, f3[1] + f3[2] + __class__.distances["F2-F3"][1])
                f2b2sel = []
                for f2_i in range(f2_s, f2_e):
                    f2 = f2b2s[f2_i]
                    if f2[1] < (f3[1] + f3[2] + __class__.distances["F2-F3"][0]):
                        continue
                    if f2[1] > (f3[1] + f3[2] + __class__.distances["F2-F3"][1]):
                        continue
                    f2b2sel += [f2]
                f2b2sel = sorted(f2b2sel, key=lambda x: max(x[0].dg5, x[0].dg3))
                for f2 in f2b2sel:
                    if f2[1] in phash['F2']:
                        continue
                    phash['F2'][f2[1]] = f2
                    f1_s = binarySearch(f1b1s, f2[1] + f2[2] + __class__.distances["loop F1c-F2"][0])
                    f1_e = binarySearch(f1b1s, f2[1] + f2[2] + __class__.distances["loop F1c-F2"][1])
                    f1b1sel = []
                    for f1_i in range(f1_s, f1_e):
                        f1 = f1b1s[f1_i]
                        if f1[1] < (f2[1] + f2[2] + __class__.distances["loop F1c-F2"][0]):
                            continue
                        if f1[1] > (f2[1] + f2[2] + __class__.distances["loop F1c-F2"][1]):
                            continue
                        f1b1sel += [f1]
                    f1b1sel = sorted(f1b1sel, key=lambda x: max(x[0].dg5, x[0].dg3))
                    for f1 in f1b1sel:
                        if f1[1] in phash['F1']:
                            continue
                        phash['F1'][f1[1]] = f1
                        b1_s = binarySearch(f1b1s, f1[1] + f1[2] + __class__.distances["F1c-B1c"][0])
                        b1_e = binarySearch(f1b1s, f1[1] + f1[2] + __class__.distances["F1c-B1c"][1])
                        b1f1sel = []
                        for b1_i in range(b1_s, b1_e):
                            b1 = f1b1s[b1_i]
                            if b1[1] < (f1[1] + f1[2] + __class__.distances["F1c-B1c"][0]):
                                continue
                            if b1[1] > (f1[1] + f1[2] + __class__.distances["F1c-B1c"][1]):
                                continue
                            b1f1sel += [b1]
                        b1f1sel = sorted(b1f1sel, key=lambda x: max(x[0].dg5, x[0].dg3))
                        for b1 in b1f1sel:
                            if b1[1] in phash['B1']:
                                continue
                            phash['B1'][b1[1]] = b1
                            b2_s = binarySearch(f2b2s, f2[1] + f2[2] + __class__.distances["F2-B2"][0])
                            b2_e = binarySearch(f2b2s, f2[1] + f2[2] + __class__.distances["F2-B2"][1])
                            b2f2sel = []
                            for b2_i in range(b2_s, b2_e):
                                b2 = f2b2s[b2_i]
                                if b2[1] < (b1[1] + b1[2] + __class__.distances["loop F1c-F2"][0]):
                                    continue
                                if b2[1] > (b1[1] + b1[2] + __class__.distances["loop F1c-F2"][1]):
                                    continue
                                if b2[1] < (f2[1] + f2[2] + __class__.distances["F2-B2"][0]):
                                    continue
                                if b2[1] > (f2[1] + f2[2] + __class__.distances["F2-B2"][1]):
                                    continue
                                b2f2sel += [b2]
                            b2f2sel = sorted(b2f2sel, key=lambda x: max(x[0].dg5, x[0].dg3))
                            for b2 in b2f2sel:
                                if b2[1] in phash['B2']:
                                    continue
                                phash['B2'][b2[1]] = b2
                                b3_s = binarySearch(f3b3s, b2[1] + b2[2] + __class__.distances["F2-F3"][0])
                                b3_e = binarySearch(f3b3s, b2[1] + b2[2] + __class__.distances["F2-F3"][1])
                                b3f3sel = []
                                for b3_i in range(b3_s, b3_e):
                                    b3 = f3b3s[b3_i]
                                    if b3[1] < (b2[1] + b2[2] + __class__.distances["F2-F3"][0]):
                                        continue
                                    if b3[1] > (b2[1] + b2[2] + __class__.distances["F2-F3"][1]):
                                        continue
                                    b3f3sel += [b3]
                                b3f3sel = sorted(b3f3sel, key=lambda x: max(x[0].dg5, x[0].dg3))
                                for b3 in b3f3sel:
                                    if b3[1] in phash['B3']:
                                        continue
                                    phash['B3'][b3[1]] = b3
                                    pset = LAMP([f3, f2, f1, b1, b2, b3])
                                    return pset
                return None
            res2 = getPrimerSet(f3)
            res3 = getPrimerSet(f3)
            for res1 in [res2, res3]:
                if res1 is None:
                    continue
                res += [res1]
                if (len(res) % 100) == 0:
                    print(len(res))
                lf_s = binarySearch(lflbs, res1.pset[1][1] + res1.pset[1][2])
                lf_e = binarySearch(lflbs, res1.pset[2][1])
                lfsel = []
                for lf_i in range(lf_s, lf_e):
                    lf = lflbs[lf_i]
                    if lf[1] < (res1.pset[1][1] + res1.pset[1][2]):
                        continue
                    if (lf[1] + lf[2]) > (res1.pset[2][1]):
                        continue
                    lfsel += [lf]
                lfsel = sorted(lfsel, key=lambda x: max(x[0].dg5, x[0].dg3))
                for lf in lfsel:
                    if lf[1] in phash['LF']:
                        continue
                    #phash['LF'][lf[1]] = lf
                    res1.addLF(lf)
                    break
                lb_s = binarySearch(lflbs, res1.pset[3][1] + res1.pset[3][2])
                lb_e = binarySearch(lflbs, res1.pset[4][1])
                lbsel = []
                for lb_i in range(lb_s, lb_e):
                    lb = lflbs[lb_i]
                    if lb[1] < (res1.pset[3][1] + res1.pset[3][2]):
                        continue
                    if (lb[1] + lb[2]) > (res1.pset[4][1]):
                        continue
                    lbsel += [lb]
                lbsel = sorted(lbsel, key=lambda x: max(x[0].dg5, x[0].dg3))
                for lb in lbsel:
                    if lb[1] in phash['LB']:
                        continue
                    #phash['LB'][lb[1]] = lb
                    res1.addLB(lb)
                    break
        return res

from dna_features_viewer import GraphicFeature, GraphicRecord
colordb = {'F3':'#f5b183', 'F3c': '#ee7c32',
           'B3':'#8599b1', 'B3c': '#a5bfdd',
           'F2':'#93d4ef', 'F2c':'#02b1f1',
           'B2':'#c45910', 'B2c':'#c7937b',
           'F1c':'#ffc002', 'F1':'#ffe597',
           'B1c':'#b86fdc', 'B1':'#6f2ea0',
           'LF':'#b2b3b6', 'LFc':'#91afc6',
           'LB':'#b2b3b6', 'LBc':'#91afc6', 'DNA':'#09b0ed'}
def getC(name, c='#cffccc'):
    if name in colordb:
        return colordb[name]
    return c
class DNAtarget:
    def __init__(self, seq, structure=None, num=1):
        self.seq = seq
        self.structure = structure
        self.annotations = []
        self.num = num
        return
    def addAnnRaw(self, start, end, name):
        self.annotations += [[name, start, end]]
        return
    def addAnn(self, seq, name):
        res = [[name, match.start(), match.end()] 
                             for match in re.finditer(str(seq), str(self.seq))]
        self.annotations += res
        return
    def removercDNAAnn(self, target):
        ann = []
        for k in self.annotations:
            if type(k[0]) != str:
                if k[0] is not target:
                    ann += [k]
            else:
                ann += [k]
        self.annotations = ann
        return
    def removercDNA(self, start, end):
        ann = []
        for k in self.annotations:
            if type(k[0]) != str:
                if k[1] <= end and start <= k[2]:
                    k[0].removercDNAAnn(self)
                else:
                    ann += [k]
            else:
                ann += [k]
        self.annotations = ann
        return
    def addrcDNA(self, rcDNA, start, end, rstart, rend):
        self.removercDNA(start, end)
        rcDNA.removercDNA(rstart, rend)
        self.annotations += [[rcDNA, start, end]]
        rcDNA.annotations += [[self, rstart, rend]]
        return
    def getNMatches(self, primer, i1):
        n1 = getMatches(self.seq[i1:(i1+len(primer))],
                        primer.reverse_complement())
        return n1
    def available(self, s1, e1):
        res = 1
        for k in self.annotations:
            if type(k[0]) != str:
                if k[1] <= e1 and s1 <= k[2]:
                    res = 0
        if self.structure is not None:
            if self.structure[s1:e1].replace('.', '') != '':
                res = 0
        return res
    @staticmethod
    def getComp(str1):
        res = str1 + 'c'
        res = res.replace('cc', '')
        return res
    def addLAMPAnn(self, p1):
        lkeys = ['F3', 'B3', 'F2', 'B2', 'F1c', 'B1c', 'LF', 'LB']
        for k in lkeys:
            if p1.__dict__[k] is not None:
                self.addAnn(p1.__dict__[k].seq, k)
                n1 = self.__class__.getComp(k)
                self.addAnn(p1.__dict__[k].seq.reverse_complement(), n1)
        return
    def buildStructure(self):
        self.structure = '.' * len(self.seq)
        start = None
        end = None
        for k in self.annotations:
            if type(k[0]) == str and k[1] == 0:
                start = k
            if type(k[0]) == str and k[2] == len(self.seq):
                end = k
        if start is not None and self.available(start[1], start[2]):
            startc = self.__class__.getComp(start[0])
            list1 = [k for k in self.annotations if k[0] == startc and k[1] != start[1]]
            list1 = sorted(list1, key=lambda k: k[1])
            if len(list1) > 0:
                self.structure = '(' * start[2] + '.' * (list1[0][1] - start[2]) \
                + ')' * start[2] + self.structure[list1[0][2]:]
        if end is not None and self.available(end[1], end[2]):
            endc = self.__class__.getComp(end[0])
            list1 = [k for k in self.annotations if k[0] == endc and k[1] != end[1]]
            list1 = sorted(list1, key=lambda k: k[2], reverse=True)
            if len(list1) > 0:
                n1 = end[2] - end[1]
                n2 = end[1] - list1[0][2]
                self.structure = self.structure[0:list1[0][1]] + \
                '(' * n1 + '.' * n2 + ')' * n1
        return
    def __repr__(self):
        return f"DNAtarget(seq={len(self.seq)}, num={self.num})"
    def plotDNA(self, width=5):
        features = []
        for k in self.annotations:
            if k[1] < k[2] and type(k[0]) == str:
                features += [GraphicFeature(start=k[1], end=k[2], strand=+1, color=getC(k[0]),
                                            label=k[0])]
            if k[1] > k[2] and type(k[0]) == str:
                features += [GraphicFeature(start=k[1], end=k[2], strand=-1, color=getC(k[0]),
                                            label=k[0])]
            if type(k[0]) != str:
                features += [GraphicFeature(start=k[1], end=k[2], strand=-1, color=getC('DNA'),
                                            label='DNA')]
        record = GraphicRecord(sequence_length=len(self.seq), features=features)
        return record.plot(figure_width=width)
def getMatches(str1, str2):
    end = len(str1)
    if end > len(str2):
        end = len(str2)
    i = 0
    while (i < end and str1[i] == str2[i]):
        i += 1
    return i
def getAmplicon(target, primer, num=15):
    pc = primer[-num:].reverse_complement()
    i1 = target.seq.find(str(pc))
    if i1 < 0:
        return None
    if target.available(i1, i1+15):
        seq = primer + target.seq[0:(i1-1)].reverse_complement()
        n1 = target.getNMatches(primer, i1)
        t1 = DNAtarget(seq, None, target.num)
        t1.addrcDNA(target, len(primer) - n1, len(seq), 0, i1 + n1 - 1 )
        if target.structure is not None:
            target.structure = None
        return t1
    return None
def getSelfAmplicon(target):
    target.buildStructure()
    res = []
    if target.structure[0] == '(':
        index = 0
        for i in range(len(target.seq)):
            if (target.structure[i] == '('):
                index += 1
            if (target.structure[i] == ')'):
                index -= 1
            if index == 0:
                break
        n1 = len(target.seq[i:])
        seq = target.seq[i:].reverse_complement() + target.seq
        structure = '(' * n1 + target.structure[0:i] + ')' * n1
        t1 = DNAtarget(seq, structure, target.num)
        res += [t1]
    if target.structure[-1] == ')':
        index = 0
        for i in range(len(target.seq) - 1, -1, -1):
            if (target.structure[i] == ')'):
                index += 1
            if (target.structure[i] == '('):
                index -= 1
            if index == 0:
                break
        n1 = len(target.seq[0: (i+1)])
        seq = target.seq + target.seq[0:(i+1)].reverse_complement()
        structure = '(' * n1 + target.structure[(i+1):] + ')' * n1
        t1 = DNAtarget(seq, structure, target.num)
        res += [t1]
    return res

def getTargetSeq(name):
    name = re.sub('hCDX2.*', 'CDX2', name)
    if name == 'Actb4':
        name = 'mActb4'
    shash = {'KCNJ15': 'AGTTCCTCCAGGTAATTCTTACTCAAACTTGTACCAACTTGTTTTTGACTGACAGTGAACAGTGAGAGAGTTTTCTTCATTTTGAGGAACCCTAAACACCTATCTTTCCCAAGGCAACCTGTCTGGACTGAGCATTTCTCTGACTTGACATAACTTCCCATCCAGCCAGGAGTCTGCACTCTTCAGTCTTTGCAGGCAGTAGCAGAATCCCATGGTAGCCAGGTGGGTGAAGGGGAGCGAGGACGTTCTACCTGCCTTGAAGAAGACACCTGACCTGCGGAGTGAGTGACCAGTGTTTCCAGAGCCTGGCAATGGATGCCATTCACATCGGCATGTCCAGCACCCCCCTGGTGAAGCACACTGCTGGGGCTGGGCTCAAGGCCAACAGACCCCGCGTCATGTCCAAGAGTGGGCACAGCAACGTGAGAATTGACAAAGTGGATGGCATATACCTACTCTACCTGCAAGACCTGTGGACCACAGTTATCGACATGAAGTGGAGATACAAACTCACCCTGTTCGCTGCCACTTTTGTGATGACCTGGTTCCTTTTTGGAGTCATCTACTATGCCATCGCGTTTATTCATGGGGACTTAGAACCCGGTGAGCCCATTTCAAATCATACCCCCTGCATCATGAAAGTGGACTCTCTCACTGGGGCGTTTCTCTTTTCCCTGGAATCCCAGACAACCATTGGCTATGGAGTCCGTTCC',
             'PTPRC': 'AAATCTGTGCTCAGTACTGGGGAGAAGGAAAGCAAACATATGGAGATATTGAAGTTGACCTGAAAGACACAGACAAATCTTCAACTTATACCCTTCGTGTCTTTGAACTGAGACATTCCAAGAGGAAAGACTCTCGAACTGTGTACCAGTACCAATATACAAACTGGAGTGTGGAGCAGCTTCCTGCAGAACCCAAGGAATTAATCTCTATGATTCAGGTCGTCAAACAAAAACTTCCCCAGAAGAATTCCTCTGAAGGGAACAAGCATCACAAGAGTACACCTCTACTCATTCACTGCAGGGATGGATCTCAGCAAACGGGAATATTTTGTGCTTTGTTAAATCTCTTAGAAAGTGCGGAAACAGAAGAGGTAGTGGATATTTTTCAAGTGGTAAAAGCTCTACGCAAAGCTAGGCCAGGCATGGTTTCCACATTCGAGCAATATCAATTCCTATATGACGTCATTGCCAGCACCTACCCTGCTCAGAATGGACAAGTAAAGAAAAACAACCATCAAGAAGATAAAATTGAATTTGATAATGAAGTGGACAAAGTAAAGCAGGATGCTAATTGTGTTAATCCACTTGGTGCCCCAGAAAAGCTCCCTGAAGCAAAGGAACAGGCTGAAGGTTCTGAACCCACGAGTGGCACTGAGGGGCCAGAACATTCTGTCAATGGTCCTGCAAGTCCAGCTTTAAATCAAGGTTCATAGGAAAAGACATAAATGAGGAAACTCCAAACCTCCTGTTAGCTGTTATTTCTATTTTTGTAGAAGTAGGAAGTGAAAATAGGTATACAGTGGATTAATTAAATGCAGCGAACCAATATTTGTAGAAGGGTTATATTTTACTACTGTGGAAAAATATTTAAGATAGTTTTGCCAGAACAGTTTGTACAGACGTATGCTTATTTTAAAATTTTATCTCTTATTCAGTAAAAAACAACTTCTTTGTAATCGTTATGTGTGTATATGTATGTGTGTATGGGTGTGTGTTTGTGTGAGAGACAGAGAAAGAGAGAGAATTCTTTCAAGTGAATCTAAAAGCTTTTGCTTTTCCTTTGTTTTTATGAAGAAAAAATACATTTTATATTAGAAGTGTTAACTTAGCTTGAAGGATCTGTTTTTAAAAATCATAAACTGTGTGCAGACTCAATAAAATCATGTACATTTCTGAAATGACCTCAAGATGTCCTCCTTGTTCTACTCATATATATCTATCTTATATAGTTTACTATTTTACTTCTAGAGATAGTACATAAAGGTGGTATGTGTGTGTATGCTACTACAAAAAAGTTGTTAACTAAATTAACATTGGGAAATCTTATATTCCATATATTAGCATTTAGTCCAATGTCTTTTTAAGCTTATTTAATTAAAAAATTTCCAGTGAGCTTATCATGCTGTCTTTACATGGGGTTTTCAATTTTGCATGCTCGATTATTCCCTGTACAATATTTAAAATTTATTGCTTGATACTTTTGACAACAAATTAGGTTTTGTACAATTGAACTTAAATAAATGTCATTAAAATAAATAAATGCAATATGTATTAATATTCATTGTATAAAAATAGAAGAATACAAACATATTTGTTAAATATTTACATATGAAATTTAATATAGCTATTTTTATGGAATTTTTCATTGATA',
             #NM_007393.5 983:1238
            'Actb': 'ACCTCTATGCCAACACAGTGCTGTCTGGTGGTACCACCATGTACCCAGGCATTGCTGACAGGATGCAGAAGGAGATTACTGCTCTGGCTCCTAGCACCATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTACTCTGTGTGGATCGGTGGCTCCATCCTGGCCTCACTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAGGAGTACGATGAGTCCGGCCCCTCCATCGTGCACCGCAAGTGCTTCTAGG',
             #NM_007393.5 837:1181
            'Actb2': 'TGACGGCCAGGTCATCACTATTGGCAACGAGCGGTTCCGATGCCCTGAGGCTCTTTTCCAGCCTTCCTTCTTGGGTATGGAATCCTGTGGCATCCATGAAACTACATTCAATTCCATCATGAAGTGTGACGTTGACATCCGTAAAGACCTCTATGCCAACACAGTGCTGTCTGGTGGTACCACCATGTACCCAGGCATTGCTGACAGGATGCAGAAGGAGATTACTGCTCTGGCTCCTAGCACCATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTACTCTGTGTGGATCGGTGGCTCCATCCTGGCCTCACTGTCCACCTTCCAGCAGATGTGGATCA',
             #NM_007393.5 837:1238
            'Actb3': 'TGACGGCCAGGTCATCACTATTGGCAACGAGCGGTTCCGATGCCCTGAGGCTCTTTTCCAGCCTTCCTTCTTGGGTATGGAATCCTGTGGCATCCATGAAACTACATTCAATTCCATCATGAAGTGTGACGTTGACATCCGTAAAGACCTCTATGCCAACACAGTGCTGTCTGGTGGTACCACCATGTACCCAGGCATTGCTGACAGGATGCAGAAGGAGATTACTGCTCTGGCTCCTAGCACCATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTACTCTGTGTGGATCGGTGGCTCCATCCTGGCCTCACTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAGGAGTACGATGAGTCCGGCCCCTCCATCGTGCACCGCAAGTGCTTCTAGG',
             
             #mm39_dna range=chr5:142890823-142891340 5'pad=0 3'pad=0 strand=- repeatMasking=none
             #1:22 exon 22:109 intron 109:349 exon 349:518 intron
             'mActb4': 'CGTGGGCCGCCCTAGGCACCAGGTAAGTGACCTGTTACTTTGGGAGTGGCAAGCCTGGGGTTTTCTTGGGGATCGATGCCGGTGCTAAGAAGGCTGTTCCCTTCCACAGGGTGTGATGGTGGGAATGGGTCAGAAGGACTCCTATGTGGGTGACGAGGCCCAGAGCAAGAGAGGTATCCTGACCCTGAAGTACCCCATTGAACATGGCATTGTTACCAACTGGGACGACATGGAGAAGATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCCCCTGAGGAGCACCCTGTGCTGCTCACCGAGGCCCCCCTGAACCCTAAGGCCAACCGTGAAAAGATGACCCAGGTCAGTATCCCGGGTAACCCTTCTCTTTGGCCAGCTTCTCAGCCACGCCCTTTCTCAATTGTCTTTCTTCTGCCGTTCTCCCATAGGACTCCCTTCTATGAGCTGAGTCTCCCTTGGATCTTTGCAGTTTCTGCTCTTTCCCAGACGAGGTCTTTTTTTCTCTCAATTG',
             #mm39_dna range=chr5:142889439-142889847 5'pad=0 3'pad=0 strand=- repeatMasking=none
             #1:26 exon 26:151 intron 151:409 exon
            'mActb5': 'TGGCTCCTAGCACCATGAAGATCAAGGTAAGCTAAGCATCCTTAGCTTGGTGAGGGTGGGCCCTGTGGTTGTCAGAGCAACCTTCTAGGTTTAAGGGGAATCCCAGCACCCAGAGAGCTCACCATTCACCATCTTGTCTTGCTTTCTTCAGATCATTGCTCCTCCTGAGCGCAAGTACTCTGTGTGGATCGGTGGCTCCATCCTGGCCTCACTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAGGAGTACGATGAGTCCGGCCCCTCCATCGTGCACCGCAAGTGCTTCTAGGCGGACTGTTACTGAGCTGCGTTTTACACCCTTTCTTTGACAAAACCTAACTTGCGCAGAAAAAAAAAAAATAAGAGACAACATTGGCATGGCTTTGTTTTTTTAAATTTTTTT',
             #NM_007673.3[713..1059]
            'Cdx2': 'CGCCGCCGCCGAACAGCTGTCCCCCAGCGGCCAGCGGCGAAACCTGTGCGAGTGGATGCGGAAGCCCGCGCAGCAGTCCCTAGGAAGCCAAGTGAAAACCAGGACAAAAGACAAATACCGGGTGGTGTACACAGACCATCAGCGGCTGGAGCTGGAGAAGGAGTTTCACTTTAGTCGATACATCACCATCAGGAGGAAAAGTGAGCTGGCTGCCACACTTGGGCTCTCCGAGAGGCAGGTTAAAATTTGGTTTCAGAACCGCAGAGCCAAGGAGAGGAAAATCAAGAAGAAGCAGCAGCAGCAACAGCAGCAGCAGCAACAACAGCCTCCACAGCCGCCG',
             #NM_007673.3[713..1059]
            'mCdx2': 'CGCCGCCGCCGAACAGCTGTCCCCCAGCGGCCAGCGGCGAAACCTGTGCGAGTGGATGCGGAAGCCCGCGCAGCAGTCCCTAGGAAGCCAAGTGAAAACCAGGACAAAAGACAAATACCGGGTGGTGTACACAGACCATCAGCGGCTGGAGCTGGAGAAGGAGTTTCACTTTAGTCGATACATCACCATCAGGAGGAAAAGTGAGCTGGCTGCCACACTTGGGCTCTCCGAGAGGCAGGTTAAAATTTGGTTTCAGAACCGCAGAGCCAAGGAGAGGAAAATCAAGAAGAAGCAGCAGCAGCAACAGCAGCAGCAGCAACAACAGCCTCCACAGCCGCCG',
             #NM_001101.5[660..901].flat
            'hACTB': 'CTCACCGAGCGCGGCTACAGCTTCACCACCACGGCCGAGCGGGAAATCGTGCGTGACATTAAGGAGAAGCTGTGCTACGTCGCCCTGGACTTCGAGCAAGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTACGAGCTGCCTGACGGCCAGGTCATCACCATTGGCAATGAGCGGTTCCGCTGCCCTGAGGCACTCTTCCAGCCTTCCTTCCTGGGCATGGAGTCCTGTG',
             #NM_001101.5[270..511].flat
            'ACTB2': 'GGCATCCTCACCCTGAAGTACCCCATCGAGCACGGCATCGTCACCAACTGGGACGACATGGAGAAAATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCTCCCGAGGAGCACCCCGTGCTGCTGACCGAGGCCCCCCTGAACCCCAAGGCCAACCGCGAGAAGATGACCCAGATCATGTTTGAGACCTTCAACACCCCAGCCATGTACGTTGCTATCCAGGCTGTGCTATCCCTGT',
             #NM_001265.6[728..1311].flat
            'CDX2': 'CAACCCCGGCCCTCCTGGGCCCGCCGCCACCGCTGCCGCCGAGCAGCTGTCTCCCGGCGGCCAGCGGCGGAACCTGTGCGAGTGGATGCGGAAGCCGGCGCAGCAGTCCCTCGGCAGCCAAGTGAAAACCAGGACGAAAGACAAATATCGAGTGGTGTACACGGACCACCAGCGGCTGGAGCTGGAGAAGGAGTTTCACTACAGTCGCTACATCACCATCCGGAGGAAAGCCGAGCTAGCCGCCACGCTGGGGCTCTCTGAGAGGCAGGTTAAAATCTGGTTTCAGAACCGCAGAGCAAAGGAGAGGAAAATCAACAAGAAGAAGTTGCAGCAGCAACAGCAGCAGCAGCCACCACAGCCGCCTCCGCCGCCACCACAGCCTCCCCAGCCTCAGCCAGGTCCTCTGAGAAGTGTCCCAGAGCCCTTGAGTCCGGTGTCTTCCCTGCAAGCCTCAGTGCCTGGCTCTGTCCCTGGGGTTCTGGGGCCAACTGGGGGGGTGCTAAACCCCACCGTCACCCAGTGACCCACCGGGTTCTGCAGCGGCAGAGCAATTCCAGGCTGAGCCATGAGGAGCGTGGACTCT',
             #mm10_dna range=chr5:142904346-142904780 5'pad=0 3'pad=0 strand=- repeatMasking=none
            'mActb-DNA': 'CATGTTTGAGACCTTCAACACCCCAGCCATGTACGTAGCCATCCAGGCTGTGCTGTCCCTGTATGCCTCTGGTCGTACCACAGGCATTGTGATGGACTCCGGAGACGGGGTCACCCACACTGTGCCCATCTACGAGGGCTATGCTCTCCCTCACGCCATCCTGCGTCTGGACCTGGCTGGCCGGGACCTGACAGACTACCTCATGAAGATCCTGACCGAGCGTGGCTACAGCTTCACCACCACAGCTGAGAGGGAAATCGTGCGTGACATCAAAGAGAAGCTGTGCTATGTTGCTCTAGACTTCGAGCAGGAGATGGCCACTGCCGCATCCTCTTCCTCCCTGGAGAAGAGCTATGAGCTGCCTGACGGCCAGGTCATCACTATTGGCAACGAGCGGTTCCGATGCCCTGAGGCTCTTTTCCAGCCTTCCTTCTT',
             #NM_003332.4 1:575
            'hTYROBP-RNA': 'ACTTGCCTGGACGCTGCGCCACATCCCACCGGCCCTTACACTGTGGTGTCCAGCAGCATCCGGCTTCATGGGGGGACTTGAACCCTGCAGCAGGCTCCTGCTCCTGCCTCTCCTGCTGGCTGTAAGTGGTCTCCGTCCTGTCCAGGCCCAGGCCCAGAGCGATTGCAGTTGCTCTACGGTGAGCCCGGGCGTGCTGGCAGGGATCGTGATGGGAGACCTGGTGCTGACAGTGCTCATTGCCCTGGCCGTGTACTTCCTGGGCCGGCTGGTCCCTCGGGGGCGAGGGGCTGCGGAGGCAGCGACCCGGAAACAGCGTATCACTGAGACCGAGTCGCCTTATCAGGAGCTCCAGGGTCAGAGGTCGGATGTCTACAGCGACCTCAACACACAGAGGCCGTATTACAAATGAGCCCGAATCATGACAGTCAGCAACATGATACCTGGATCCAGCCATTCCTGAAGCCCACCCTGCACCTCATTCCAACTCCTACCGCGATACAGACCCACAGAGTGCCATCCCTGAGAGACCAGACCGCTCCCCAATACTCTCCTAAAATAAACATGAAGCACAAAAA'}
    if name in shash:
        return shash[name]
    return None

def getPrimersAll(cfile, num=93):
    dfp = pd.read_excel(cfile, engine='openpyxl')
    list3 = []
    name, id1, target1 = None, None, None
    f3, f2, f1, b1, b2, b3 = [None] * 6
    lf, lb = None, None
    for i in dfp.index:
        if not pd.isna(dfp['Name'][i]):
            if pd.isna(dfp['Sequence'][i]):
                name = dfp['Name'][i]
                continue
            if dfp['Name'][i].endswith('F3'):
                id1 = re.sub("[_-]*F3", "", dfp['Name'][i])
                target1 = dfp['Target'][i]
                pos = 0
                f3 = [Seq(dfp['Sequence'][i]), pos, len(dfp['Sequence'][i])]
                f3[0].id = dfp['Name'][i]
            if dfp['Name'][i].endswith('B3'):
                pos = 0
                b3 = [Seq(dfp['Sequence'][i]).reverse_complement(), pos, len(dfp['Sequence'][i])]
                b3[0].id = dfp['Name'][i]
            if dfp['Name'][i].endswith('LF'):
                pos = 0
                lf = [Seq(dfp['Sequence'][i]).reverse_complement(), pos, len(dfp['Sequence'][i])]
                lf[0].id = dfp['Name'][i]
            if dfp['Name'][i].endswith('LB'):
                pos = 0
                lb = [Seq(dfp['Sequence'][i]), pos, len(dfp['Sequence'][i])]
                lb[0].id = dfp['Name'][i]
            if dfp['Name'][i].endswith('FIP'):
                f1c, f2s = splitPrimers(target1, dfp['Sequence'][i])
                pos = 0
                f1 = [Seq(f1c).reverse_complement(), pos, len(f1c)]
                f2 = [Seq(f2s), pos, len(f2s)]
                f1[0].id = dfp['Name'][i]
                f2[0].id = dfp['Name'][i]
            if dfp['Name'][i].endswith('BIP'):
                b1c, b2s = splitPrimers(Seq(target1).reverse_complement(),
                                             dfp['Sequence'][i])
                pos = 0
                b1 = [Seq(b1c), pos, len(b1c)]
                b2 = [Seq(b2s).reverse_complement(), pos, len(b2s)]
                b1[0].id = dfp['Name'][i]
                b2[0].id = dfp['Name'][i]
        if pd.isna(dfp['Name'][i]) and f3 is None:
            name, id1, target1 = None, None, None
            f3, f2, f1, b1, b2, b3 = [None] * 6
            lf, lb = None, None
        if pd.isna(dfp['Name'][i]) and f3 is not None:
            pset = LAMP([f3, f2, f1, b1, b2, b3])
            pset.F3.id = f3[0].id
            pset.F2.id = f2[0].id
            pset.F1c.id = f1[0].id
            pset.B1c.id = b1[0].id
            pset.B2.id = b2[0].id
            pset.B3.id = b3[0].id
            if lf is not None:
                pset.addLF(lf)
                pset.LF.id = lf[0].id
                pset.lfdata = lf
            if lb is not None:
                pset.addLB(lb)
                pset.LB.id = lb[0].id
                pset.lbdata = lb
            pset.name = name
            pset.id1 = id1
            pset.target1 = target1
            list3 += [pset]
            name, id1, target1 = None, None, None
            f3, f2, f1, b1, b2, b3 = [None] * 6
            lf, lb = None, None
        if (len(list3) == num):
            break
    return list3

def getPrimersMatching(list3, seq, mset):
    minlen = LAMP.distances["F2-B2"][0] + 2 * LAMP.lengths["F3/B3"][0]
    m1 = np.min(np.array([LAMP.lengths[k] for k in LAMP.lengths]))
    m2 = np.max(np.array([LAMP.lengths[k] for k in LAMP.lengths]))
    n = len(seq)
    primers = []
    f3b3s = []
    lflbs = []
    fips = []
    bips = []
    phash = {'F3':{},'LF':{},'All':{}, 'FIP':{}, 'BIP':{}}
    for l1 in list3:
        lkeys = ['F3', 'B3', 'LF', 'LB']
        chash = {'F3':0, 'B3':1, 'LF':1, 'LB':0}
        for k in lkeys:
            if l1.__dict__[k] is None:
                continue
            if chash[k] == 0:
                primer = LAMP.getSeqParamaters(l1.__dict__[k].seq)
                i1 = seq.find(primer)
            else:
                primer = LAMP.getSeqParamaters(l1.__dict__[k].seq.reverse_complement())
                i1 = seq.find(primer)
            if i1 < 0:
                continue
            primer.id = l1.__dict__[k].id + '-' + str(l1.id)
            j = len(primer)
            pos = i1
            if LAMP.isStable(primer) and pos not in phash['All']:
                primers += [[primer, pos, j]]
                phash['All'][pos] = [primer, pos, j]
            if LAMP.isValid(primer, "F3/B3") and pos not in phash['F3']:
                f3b3s += [[primer, pos, j]]
                phash['F3'][pos] = [primer, pos, j]
            if LAMP.isValid(primer, "LF/LB") and pos not in phash['LF']:
                lflbs += [[primer, pos, j]]
                phash['LF'][pos] = [primer, pos, j]
        pos1 = getPos(seq, l1.__dict__['F1c'].seq)
        pos2 = getPos(seq, l1.__dict__['F2'].seq)
        if pos1 >= 0 and pos2 >= 0:
            pos = f"{pos1}-{pos2}"
            if pos not in phash['FIP']:
                f1 = l1.__dict__['F1c'].seq.reverse_complement()
                f2 = l1.__dict__['F2'].seq
                f1.id = l1.__dict__['F1c'].id + '-' + str(l1.id)
                f2.id = l1.__dict__['F2'].id + '-' + str(l1.id)
                fips += [[f1, pos1, len(f1), f2, pos2, len(f2)]]
                phash['FIP'][pos] = 1
        pos1 = getPos(seq, l1.__dict__['B1c'].seq)
        pos2 = getPos(seq, l1.__dict__['B2'].seq)
        if pos1 >= 0 and pos2 >= 0:
            pos = f"{pos1}-{pos2}"
            if pos not in phash['BIP']:
                b1 = l1.__dict__['B1c'].seq
                b2 = l1.__dict__['B2'].seq.reverse_complement()
                b1.id = l1.__dict__['B1c'].id + '-' + str(l1.id)
                b2.id = l1.__dict__['B2'].id + '-' + str(l1.id)
                bips += [[b1, pos1, len(b1), b2, pos2, len(b2)]]
                phash['BIP'][pos] = 1
    list4 = []
    primers = sorted(primers, key=lambda x: x[1])
    f3b3s = sorted(f3b3s, key=lambda x: x[1])
    lflbs = sorted(lflbs, key=lambda x: x[1])
    pdata = [primers, f3b3s, lflbs, fips, bips]
    print(len(primers), len(f3b3s), len(lflbs), len(fips), len(bips))
    for f3 in f3b3s:
        if mset.F3 is not None and mset.F3.seq != f3[0]:
            continue
        if f3[1] > (n - minlen):
            continue
        for fip in fips:
            f1 = fip[:3]
            f2 = fip[3:]
            if mset.F1c is not None and mset.F1c.seq != f1[0].reverse_complement():
                continue
            if mset.F2 is not None and mset.F2.seq != f2[0]:
                continue
            if f2[1] < (f3[1] + f3[2] + LAMP.distances["F2-F3"][0]):
                continue
            if f2[1] > (f3[1] + f3[2] + LAMP.distances["F2-F3"][1]):
                continue
            for bip in bips:
                b1 = bip[:3]
                b2 = bip[3:]
                if mset.B1c is not None and mset.B1c.seq != b1[0]:
                    continue
                if mset.B2 is not None and mset.B2.seq != b2[0].reverse_complement():
                    continue
                if b1[1] < (f1[1] + f1[2] + LAMP.distances["F1c-B1c"][0]):
                    continue
                if b1[1] > (f1[1] + f1[2] + LAMP.distances["F1c-B1c"][1]):
                    continue
                if b2[1] < (f2[1] + f2[2] + LAMP.distances["F2-B2"][0]):
                    continue
                if b2[1] > (f2[1] + f2[2] + LAMP.distances["F2-B2"][1]):
                    continue
                for b3 in f3b3s:
                    if mset.B3 is not None and mset.B3.seq != b3[0].reverse_complement():
                        continue
                    if b3[1] < (b2[1] + b2[2] + LAMP.distances["F2-F3"][0]):
                        continue
                    if b3[1] > (b2[1] + b2[2] + LAMP.distances["F2-F3"][1]):
                        continue
                    pset = LAMP([f3, f2, f1, b1, b2, b3])
                    pset.F3.id = f3[0].id
                    pset.F2.id = f2[0].id
                    pset.F1c.id = f1[0].id
                    pset.B1c.id = b1[0].id
                    pset.B2.id = b2[0].id
                    pset.B3.id = b3[0].id
                    list4 += [pset]
                #break
            #break
        #break
    for res1 in list4:
        lfsel = []
        for lf in lflbs:
            if mset.LF is not None and mset.LF.seq != lf[0].reverse_complement():
                continue
            if lf[1] < (res1.pset[1][1] + res1.pset[1][2]):
                continue
            if (lf[1] + lf[2]) > (res1.pset[2][1]):
                continue
            lfsel += [lf]
        lbsel = []
        for lb in lflbs:
            if mset.LB is not None and mset.LB.seq != lb[0]:
                continue
            if lb[1] < (res1.pset[3][1] + res1.pset[3][2]):
                continue
            if (lb[1] + lb[2]) > (res1.pset[4][1]):
                continue    
            lbsel += [lb]
        if len(lfsel) > 0:
            res1.addLF(lfsel[0])
        if len(lbsel) > 0:
            res1.addLB(lbsel[0])
        res1.lfdata = lfsel
        res1.lbdata = lbsel
        res1.pdata = pdata
    return list4