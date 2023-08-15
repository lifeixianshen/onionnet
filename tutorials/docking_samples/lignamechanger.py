#!/usr/bin/env python

import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict


class rewritePDB(object):
    """
    Modify pdb file by changing atom indexing, resname, res sequence number and chain id

    Parameters
    ----------

    Attributes
    ----------

    """
    def __init__(self, inpdb):
        self.pdb = inpdb

    def pdbRewrite(self, input, output, chain, atomStartNdx, resStartNdx):
        """
        change atom id, residue id and chain id
        :param input: str, input pdb file
        :param output: str, output file name
        :param chain: str, chain id
        :param atomStartNdx: int,
        :param resStartNdx: int
        :return:
        """
        resseq = int(resStartNdx)
        atomseq = int(atomStartNdx)
        chainname = chain

        with open(output,'w') as newfile:
            resseq_list = []

            try:
                with open(input) as lines:
                    for s in lines:
                        if "ATOM" in s and len(s.split()) > 6:
                            atomseq += 1
                            newline = s
                            newline = self.atomSeqChanger(newline, atomseq)
                            newline = self.chainIDChanger(newline, chainname)
                            if resseq_list:
                                if resseq_list[-1] != int(s[22:26].strip()):
                                    resseq += 1
                            newline = self.resSeqChanger(newline, resseq)
                            resseq_list.append(int(s[22:26].strip()))
                            newfile.write(newline)
                        else:
                            newfile.write(s)
            except FileExistsError:
                print(f"File {self.pdb} not exist")

        return 1

    def resSeqChanger(self, inline, resseq):
        resseqstring = " "*(4 - len(str(resseq)))+str(resseq)
        return inline[:22] + resseqstring + inline[26:]

    def atomSeqChanger(self, inline, atomseq):
        atomseqstring = " " * (5 - len(str(atomseq))) + str(atomseq)
        return inline[:6] + atomseqstring + inline[11:]

    def resNameChanger(self, inline, resname):
        resnamestr = " " * (4 - len(str(resname))) + str(resname)
        return inline[:16] + resnamestr + inline[20:]

    def chainIDChanger(self, inline, chainid):
        return inline[:21] + str(chainid) + inline[22:]

    def atomNameChanger(self, inline, new_atom_name):
        return inline[:12] + "%4s" % new_atom_name + inline[16:]

    def combinePDBFromLines(self, output, lines):
        """
        combine a list of lines to a pdb file

        Parameters
        ----------
        output
        lines

        Returns
        -------

        """

        with open(output, "wb") as tofile :
            tmp = map(lambda x: tofile.write(x), lines)

        return 1

    def swampPDB(self, input, atomseq_pdb, out_pdb, chain="B"):
        """
        given a pdb file (with coordinates in a protein pocket), but with wrong atom
        sequence order, try to re-order the pdb for amber topology building

        Parameters
        ----------
        input:str,
            the pdb file with the correct coordinates
        atomseq_pdb:str,
            the pdb file with correct atom sequences
        out_pdb: str,
            output pdb file name
        chain: str, default is B
            the chain identifier of a molecule

        Returns
        -------

        """

        with open("temp.pdb", 'w') as tofile:
            crd_list = {}

            ln_target, ln_source = 0, 0
            # generate a dict { atomname: pdbline}
            with open(input) as lines :
                for s in [x for x in lines if ("ATOM" in x or "HETATM" in x)]:
                    crd_list[s.split()[2]] = s
                    ln_source += 1

            # reorder the crd_pdb pdblines, according to the atomseq_pdb lines
            with open(atomseq_pdb) as lines:
                for s in [x for x in lines if ("ATOM" in x or "HETATM" in x)]:
                    newline = crd_list[s.split()[2]]
                    tofile.write(newline)
                    ln_target += 1

        if ln_source != ln_target:
            print(
                f"Error: Number of lines in source and target pdb files are not equal. ({input} {atomseq_pdb})"
            )

        # re-sequence the atom index
        self.pdbRewrite(input="temp.pdb", atomStartNdx=1, chain=chain, output=out_pdb, resStartNdx=1)

        os.remove("temp.pdb")

        return None

def lig_name_change(lig_in, lig_out, lig_code):

    pio = rewritePDB(lig_in)
    with open(lig_out, "w") as tofile:
        with open(lig_in) as lines:
            for s in lines:
                if len(s.split()) and s.split()[0] in ['ATOM', 'HETATM']:
                    nl = pio.resNameChanger(s, lig_code)
                    #n2 = pio.chainIDChanger(nl, "Z")
                    tofile.write(nl)

    return None

def main():
    if len(sys.argv) <= 3:
        print("""Usage: \npython ligandnamechanger.py old_ligand.pdb new_ligand.pdb""")

    lig = sys.argv[1]
    out = sys.argv[2]

    with open(lig) as lines:
        lines = [x for x in lines if "LIG" in x]
        if not len(lines):
            lig_name_change(lig, out, "LIG")

        else:
            os.system(f"cp {lig} temp")

main()

