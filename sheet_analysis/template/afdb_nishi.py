import foldcomp
import pymol
from pymol import stored
import numpy as np
import re
import csv
import sys
args = sys.argv


import subprocess
import foldcomp
#from pymol import cmd
#import numpy as np


import re
def parse_nishi(lines):
    #print(lines)
    Helix_SSE_IDs = list()
    Helix_INIs = list()
    Helix_ENDs = list()

    Member_Sheet = list()
    Member_Description = list()

    NomC_Sheet = list()
    NomC_Description = list()

    NomR_Sheet = list()
    NomR_Description = list()

    ResPair_ResNum1 = list()
    ResPair_ResNum2 = list()
    ResPair_PorA = list()
    ResPair_PairType = list()
    ResPair_ForB = list()

    SHEETINFO_SheetID = list()
    SHEETINFO_Nstrands = list()
    SHEETINFO_Cycle = list()
    SHEETINFO_Undirected = list()
    SHEETINFO_WithBranch = list()
    SHEETINFO_Consecutive = list()
    SHEETINFO_AllPara = list()
    SHEETINFO_AllAnti = list()

    StrandPair_B1 = list()
    StrandPair_B2 = list()
    StrandPair_B1_B2 = list()
    StrandPair_Sheet = list()
    StrandPair_Dir = list()
    StrandPair_PorA = list()
    StrandPair_Jump = list()
    StrandPair_D1 = list()
    StrandPair_D2 = list()
    StrandPair_Bridge = list()
    StrandPair_Score = list()
    StrandPair_SSEsLBTS = list()
    StrandPair_NumResLBTS = list()

    SUBSTRAND_SubstrandID = list()
    SUBSTRAND_SheetID = list()
    SUBSTRAND_Ini = list()
    SUBSTRAND_End = list()

    for line in lines:
        line = re.sub('\s+', ' ', line)
        key = line.split(" ")[0]
        #print(key)
        if (key == "HELIX"):
            # SSE_ID Ini End
            SSE_ID = line.split(" ")[1]
            Ini = line.split(" ")[2]
            End = line.split(" ")[3]

            Helix_SSE_IDs.append(SSE_ID)
            Helix_INIs.append(Ini)
            Helix_ENDs.append(End)

        # Topology
        if (key == "MEMBER"):
            Sheet = line.split(" ")[1]
            Description = line.split(" ")[2]

            Member_Sheet.append(Sheet)
            Member_Description.append(Description)

        if (key == "NOMENCLATURE_C"):
            Sheet = line.split(" ")[1]
            Description = line.split(" ")[2]

            NomC_Sheet.append(Sheet)
            NomC_Description.append(Description)

        if (key == "NOMENCLATURE_R"):
            Sheet = line.split(" ")[1]
            Description = line.split(" ")[2]

            NomR_Sheet.append(Sheet)
            NomR_Description.append(Description)

        # Topology

#        if (key == "REMARK"):
#           print("remark")

        if (key == "RESIDUE_PAIR"):
            # ResNum1
            # ResNum2
            # PorA
            # Pair - type
            # ForB
            ResNum1 = line.split(" ")[1]
            ResNum2 = line.split(" ")[2]
            PorA = line.split(" ")[3]
            Pair_type = line.split(" ")[4]
            ForB = line.split(" ")[5]

            ResPair_ResNum1.append(ResNum1)
            ResPair_ResNum2.append(ResNum2)
            ResPair_PorA.append(PorA)
            ResPair_PairType.append(Pair_type)
            ResPair_ForB.append(ForB)

        if (key == "SHEET_INFO"):
            # REMARK
            # Sheet_ID
            # N_strands
            # Cycle
            # Undirected
            # With_branch
            # Consecutive
            # All_para
            # All_anti
            Sheet_ID = line.split(" ")[1]
            N_strands = line.split(" ")[2]
            Cycle = line.split(" ")[3]
            Undirected = line.split(" ")[4]
            With_branch = line.split(" ")[5]
            Consecutive = line.split(" ")[6]
            All_para = line.split(" ")[7]
            All_anti = line.split(" ")[8]

            SHEETINFO_SheetID.append(Sheet_ID)
            SHEETINFO_Nstrands.append(N_strands)
            SHEETINFO_Cycle.append(Cycle)
            SHEETINFO_Undirected.append(Undirected)
            SHEETINFO_WithBranch.append(With_branch)
            SHEETINFO_Consecutive.append(Consecutive)
            SHEETINFO_AllPara.append(All_para)
            SHEETINFO_AllAnti.append(All_anti)

        if (key == "STRAND_PAIR"):
            # REMARK
            # B1
            # B2
            # Sheet
            # Dir
            # PorA
            # Jump
            # D1
            # D2
            # Bridge
            # Score
            # SSEs_LBTS
            # NumRes_LBTS

            B1 = line.split(" ")[1]
            B2 = line.split(" ")[2]
            Sheet = line.split(" ")[3]
            Dir = line.split(" ")[4]
            PorA = line.split(" ")[5]
            Jump = line.split(" ")[6]
            D1 = line.split(" ")[7]
            D2 = line.split(" ")[8]
            Bridge = line.split(" ")[9]
            Score = line.split(" ")[10]
            SSEs_LBTS = line.split(" ")[11]
            NumRes_LBTS = line.split(" ")[12]

            StrandPair_B1.append(B1)
            StrandPair_B2.append(B2)
            B1B2 = B1 + "+" + B2
            StrandPair_B1_B2.append(B1B2)
            StrandPair_Sheet.append(Sheet)
            StrandPair_Dir.append(Dir)
            StrandPair_PorA.append(PorA)
            StrandPair_Jump.append(Jump)
            StrandPair_D1.append(D1)
            StrandPair_D2.append(D2)
            StrandPair_Bridge.append(Bridge)
            StrandPair_Score.append(Score)
            StrandPair_SSEsLBTS.append(SSEs_LBTS)
            StrandPair_NumResLBTS.append(NumRes_LBTS)

        if (key == "SUBSTRAND"):
            # REMARK
            # SubStrand_ID
            # Sheet_ID
            # Ini
            # End
            SubStrand_ID = line.split(" ")[1]
            Sheet_ID = line.split(" ")[2]
            Ini = line.split(" ")[3]
            End = line.split(" ")[4]

            SUBSTRAND_SubstrandID.append(SubStrand_ID)
            SUBSTRAND_SheetID.append(Sheet_ID)
            SUBSTRAND_Ini.append(Ini)
            SUBSTRAND_End.append(End)
        # NOMENCLATURE_R
        # REMARK
        # RESIDUE_PAIR
        # SHEET_INFO
        # STRAND_PAIR
        # SUBSTRAND

    #print(SUBSTRAND_SheetID)
    HELIX = {"SSEIDs": Helix_SSE_IDs,
             "Ini": Helix_INIs,
             "End": Helix_ENDs
             }

    MEMBER = {"Sheet": Member_Sheet,
              "Description": Member_Description
              }

    NomC = {"Sheet": NomC_Sheet,
            "Description": NomC_Description
            }

    NomR = {"Sheet": NomR_Sheet,
            "Description": NomR_Description
            }

    TOPOLOGY = {"MEMBER": MEMBER,
                "NomC": NomC,
                "NomR": NomR}

    ResPair = {"ResNum1": ResPair_ResNum1,
               "ResNum2": ResPair_ResNum2,
               "PorA": ResPair_PorA,
               "PairType": ResPair_PairType,
               "ForB": ResPair_ForB,
               }

    SHEETINFO = {"SheetID": SHEETINFO_SheetID,
                 "Nstrands": SHEETINFO_Nstrands,
                 "Cycle": SHEETINFO_Cycle,
                 "Undirected": SHEETINFO_Undirected,
                 "WithBranch": SHEETINFO_WithBranch,
                 "Consecutive": SHEETINFO_Consecutive,
                 "AllPara": SHEETINFO_AllPara,
                 "AllAnti": SHEETINFO_AllAnti,
                 }

    STRANDPAIR = {"B1": StrandPair_B1,
                  "B2": StrandPair_B2,
                  "B1B2": StrandPair_B1_B2,
                  "Sheet": StrandPair_Sheet,
                  "Dir": StrandPair_Dir,
                  "PorA": StrandPair_PorA,
                  "Jump": StrandPair_Jump,
                  "D1": StrandPair_D1,
                  "D2": StrandPair_D2,
                  "Bridge": StrandPair_Bridge,
                  "Score": StrandPair_Score,
                  "SSEsLBTS": StrandPair_SSEsLBTS,
                  "NumResLBTS": StrandPair_NumResLBTS
                  }

    SUBSTRAND = {"SubstrandID": SUBSTRAND_SubstrandID,
                 "SheetID": SUBSTRAND_SheetID,
                 "Ini": SUBSTRAND_Ini,
                 "End": SUBSTRAND_End
                 }

    NISHI = {"SubStrand": SUBSTRAND,
             "StrandPair": STRANDPAIR,
             "SheetInfo": SHEETINFO,
             "ResPair": ResPair,
             "Topology": TOPOLOGY,
             "Helix": HELIX,
             }

    return NISHI

def search_AfterWalker(nishi=None,init=None):
    # find beta strand just before the ploop
    target_i = -1
    for i in range(0,len(nishi["SubStrand"]["Ini"])):
        if (int(nishi["SubStrand"]["Ini"][i]) <= init and int(nishi["SubStrand"]["End"][i]) >= init):
            #print("HIT SubID=", nishi["SubStrand"]["SubstrandID"][i])
            target_i = i
            strandID_1 = nishi["SubStrand"]["SubstrandID"][i]
            #print("HIT SheID=",nishi["SubStrand"]["SheetID"][i])
            sheetID_1 = nishi["SubStrand"]["SheetID"][i]

    # find beta strand after the ploop and in the same sheet as Sheet1
    if (target_i>=0):
        strandID_2 = -99
        for j in range(target_i+1, len(nishi["SubStrand"]["Ini"])):
            #print(nishi["SubStrand"])
            if (nishi["SubStrand"]["SheetID"][j]==sheetID_1):
                target_j = j
                strandID_2 = nishi["SubStrand"]["SubstrandID"][j]
                sheetID_2 = nishi["SubStrand"]["SheetID"][j]
                #print("HIT SubID2=", nishi["SubStrand"]["SubstrandID"][j])
                #print("HIT SubID2=", nishi["SubStrand"]["SheetID"][j])
                break
        if (strandID_2 != -99):
            targetpair = strandID_1 + "+" + strandID_2
            for k in range(0, len(nishi["StrandPair"]["B1B2"])):
                if (nishi["StrandPair"]["B1B2"][k]==targetpair):
                    target_k = k
                    break

#    SHEETINFO = {"SheetID": SHEETINFO_SheetID,
#                 "Nstrands": SHEETINFO_Nstrands,
#                 "Cycle": SHEETINFO_Cycle,
#                 "Undirected": SHEETINFO_Undirected,
#                 "WithBranch": SHEETINFO_WithBranch,
#                 "Consecutive": SHEETINFO_Consecutive,
#                 "AllPara": SHEETINFO_AllPara,
#                 "AllAnti": SHEETINFO_AllAnti,
#                 }


            for l in range(0,len(nishi["SheetInfo"]["SheetID"])):
                if (nishi["SheetInfo"]["SheetID"][l]==sheetID_1):
                    all_para = nishi["SheetInfo"]["AllPara"][l]
                    all_anti = nishi["SheetInfo"]["AllAnti"][l]
                    cycle = nishi["SheetInfo"]["Cycle"][l]
                    undire = nishi["SheetInfo"]["Undirected"][l]
                    withbranch = nishi["SheetInfo"]["WithBranch"][l]
                    consecutive = nishi["SheetInfo"]["Consecutive"][l]
                    nstrands = nishi["SheetInfo"]["Nstrands"][l]
                    break

            return nishi["StrandPair"]["Jump"][target_k],\
                nishi["StrandPair"]["PorA"][target_k], \
                nishi["StrandPair"]["Dir"][target_k],\
                sheetID_1, \
                all_para, \
                all_anti,\
                cycle,\
                undire,\
                withbranch,\
                consecutive,\
                nstrands
        else:
            print("Not paired to latter strand")
            return -2,"neither","neither",-2,"neither","neither","neigher","neither","neither","neither","neither"
    else:
        print("Neither")
        return -1,"neither","neither",-1,"neither","neither","neigher","neither","neither","neither","neither"
        
def ppos2abegos(ppos,divideB):
    abegos=[]
    for resi in range(0,len(ppos[:,0])):
        phi   = float(ppos[resi, 0])
        psi   = float(ppos[resi, 1])
        omega = float(ppos[resi, 2])
        abegos.append(dihd2abego(phi,psi,omega,divideB=divideB))
    if (len(ppos[:,0] == len(abegos))):
        return abegos
    else:
        return None

def dihd2abego(phi=-60.0,psi=-45.0,omega=180.0,
               cisbin=30.0,A_upper=50.0,A_lower=-75.0,
               G_upper=100,G_lower=-100.0,
                P_thre = -90, divideB=False):
    if (omega <= cisbin and omega >= -cisbin):
        return "O"

    if (phi <= 0):
        if (A_lower <= psi and psi <= A_upper):
            return "A"
        else:
            if (divideB):
                if ( phi > P_thre):
                    return "P"
                else:
                    return "S"
            else:
                return "B"
    else:
        if (G_lower <= psi and psi <= G_upper):
            return "G"
        else:
            return "E"
####################################

def pdb2abego(pdb):
    pymol.cmd.delete("all")
    pymol.cmd.read_pdbstr(pdb, oname="strc")
    fasta=pymol.cmd.get_fastastr("strc")
    myfasta=fasta[8:].replace("\n","") # hard coded for 4-letter name : strc_A\n
    phipsi = pymol.cmd.phi_psi("strc")
    keys = list(phipsi)
    stored.bfactor = []
    pymol.cmd.iterate("name ca", "stored.bfactor.append(b)")
    ppos = []
    for key in keys:
        phi = phipsi[key][0]
        psi = phipsi[key][1]
        omega = 180.0  # temporary
        ppos.append([phi, psi, omega])
    ppos = np.array(ppos)
    myabego = ppos2abegos(ppos, divideB=False)
    return "".join(myabego), myfasta, stored.bfactor

def iterative_search(name,abego,fasta,bfactor,query="BBBEBBGAGAAAAA"):
    fasta_tmp = fasta[1:]
    abego_tmp = abego
    bfactor_tmp = bfactor[1:]

    seqs=list()
    inits = list()
    ends=list()
    mynames = list()
    queries = list()
    bfactors = list()
    meanbfactors = list()
    GKTs = list()
    IDs = list()
    rmsds = list()
    fastas = list()
    cumm=0

    result = re.search(query, abego_tmp)
    if (result==None):
        #print("HOGEHOGEHOGEHOGE")
        return mynames, queries, seqs, inits, ends, bfactors,meanbfactors,GKTs,IDs,rmsds,fastas
    else:
        #print("PIYO")
        # save fcz structrure
        with open("./structures/" + name + ".fcz", "bw") as fcz_file:
            cmp = foldcomp.compress(pdb_content=pdb, name=name)
            fcz_file.write(cmp)
        # save fasta

        pymol.cmd.delete("all")
        pymol.cmd.load("1BMF_A_166-179.pdb")
        pymol.cmd.read_pdbstr(pdb, oname=name)
        pymol.cmd.save("./fastas/" + name + ".fasta", name)

        while (len(fasta_tmp) > 0) :
            result = re.search(query, abego_tmp)
            if (result == None):
                break
            else:
                i = result.span()[0]
                e = result.span()[1]
                abs_i = i + 2 + cumm
                abs_e = e + 1 + cumm
                seqs.append(fasta_tmp[i:e])
                GKTs.append(fasta_tmp[i:e][8:11])
                #EBBGAGAA
                #BBBEBBGAGAAA
                #012345678
                bfactors.append(sum(bfactor_tmp[i:e]))
                meanbfactors.append(sum(bfactor) / len(bfactor))

                fastas.append(fasta)

                inits.append(abs_i)
                ends.append(abs_e)
                mynames.append(name)
                IDs.append(name.split(" ")[-1].replace("(","").replace(")",""))
                queries.append(query)
                selectname=name+"_"+str(abs_i)+"_"+str(abs_e)
                pymol.cmd.create(selectname, name + " and resi " + str(abs_i)+"-"+str(abs_e))

                rmsd=pymol.cmd.pair_fit(selectname +" and name ca","1BMF_A_166-179 and name ca")
                rmsds.append(rmsd)
                
                pymol.cmd.save("./fragments/"+name+"_"+str(abs_i)+"_"+str(abs_e)+".pdb",selectname)
                pymol.cmd.save("./fragment_fastas/"+name+"_"+str(abs_i)+"_"+str(abs_e)+".fasta",selectname)
                
                pymol.cmd.delete(selectname)
                
                fasta_tmp = fasta_tmp[e:]
                abego_tmp = abego_tmp[e:]
                bfactor_tmp = bfactor_tmp[e:]
                cumm += e

    return mynames, queries, seqs, inits, ends, bfactors, GKTs, meanbfactors, IDs, rmsds, fastas



#    return nishi["StrandPair"]["Jump"][target_k],\
#                nishi["StrandPair"]["PorA"][target_k], \
#                nishi["StrandPair"]["Dir"][target_k],\
#                sheetID_1, \
#                all_para, \
#                all_anti,
#                cycle,
#                undire,
#                withbranch,
#                consecutive,
#                nstrands


execute = "/lab/home/ksakuma/mybin/beta_sheet_analysis/bin/PROG_NAME -w target.pdb"

#option1 = "-w"
nishi_prog = [execute]
with foldcomp.open("./splitdb/uniprotv4walkerlike") as db:
  # Iterate through database
  names = list()
  queries = list()
  seqs = list()
  inits = list()
  ends = list()
  bfactors = list()
  gkts = list()
  meanbfactors = list()
  IDs = list()
  rmsds = list()
  fastas = list()
  jumps = list()
  PorAs = list()
  dires = list()
  SheetIDs = list()
  all_paras = list()
  all_antis = list()
  cycles = list()
  undires = list()
  withbranchs = list()
  consecutives= list()
  nstrandss=list()


  count = 0
  pcount =0
  import time
  for (name_tmp, pdb) in db:
      name=name_tmp.replace(",","").replace("/","").replace(" ","_")
      abego,fasta,b=pdb2abego(pdb)
      ####### sheet analysis #########
      #print(pdb)
      pdb_file=open("./target.pdb",mode="w")
      pdb_file.write(pdb)
      pdb_file.close()
      proc = subprocess.run(nishi_prog, stdout=subprocess.PIPE,shell=True)
      lines = proc.stdout.decode("utf8").split("\n")
      print(name)
      nishi_result = parse_nishi(lines)
     #################################

      count += 1
      if (count % 1000==0):
          print(name)
          print(count)
          print(len(IDs))

      if (len(IDs) >= 1000):
          pcount=pcount+1
          f = open('afdb_ploop_'+str(pcount).zfill(10)+".csv", 'w')
          data = list()
          data.append(names)
          data.append(queries)
          data.append(seqs)
          data.append(inits)
          data.append(ends)
          data.append(bfactors)
          data.append(gkts)
          data.append(meanbfactors)
          data.append(IDs)
          data.append(rmsds)
          data.append(fastas)

          data.append(jumps)
          data.append(PorAs)
          data.append(dires)
          data.append(SheetIDs)
          data.append(all_paras)
          data.append(all_antis)

          data.append(cycles)
          data.append(undires)
          data.append(withbranchs)
          data.append(consecutives)
          data.append(nstrandss)


          data2 = list(zip(*data))
          writer = csv.writer(f)
          writer.writerows(data2)
          f.close()

          names = list()
          queries = list()
          seqs = list()
          inits = list()
          ends = list()
          bfactors = list()
          gkts=list()
          meanbfactors=list()
          IDs = list()
          rmsds = list()
          fastas = list()
          jumps= list()
          PorAs = list()
          dires = list()
          SheetIDs = list()
          all_paras = list()
          all_antis = list()
          cycles = list()
          undires = list()
          withbranchs = list()
          consecutives= list()
          nstrandss=list()

      for qabego in ["BBBBBBGAGAAAAA","BBBEBBGAGAAAAA","BBBABBGAGAAAAA","BBBGBBGAGAAAAA"]:
          myname, query, seq, init, end, bfactor, gkt, meanbfactor, ID, rmsd, fst = iterative_search(name=name,fasta=fasta,abego=abego,query=qabego,bfactor=b)

#          print("init",init)
#          print("end",end)
          for ini in init:
              jump,PorA,dire,SheetID,all_para,all_anti,cycle,undire,withbranch,consecutive,nstrands=search_AfterWalker(nishi=nishi_result,init=ini)
              jumps.append(jump)
              dires.append(dire)
              PorAs.append(PorA)
              SheetIDs.append(SheetID)
              all_paras.append(all_para)
              all_antis.append(all_anti)
              cycles.append(cycle)
              undires.append(undire)
              withbranchs.append(withbranch)
              consecutives.append(consecutive)
              nstrandss.append(nstrands)

          names.extend(myname)
          queries.extend(query)
          seqs.extend(seq)
          inits.extend(init)
          ends.extend(end)
          bfactors.extend(bfactor)
          gkts.extend(gkt)
          meanbfactors.extend(meanbfactor)
          IDs.extend(ID)
          rmsds.extend(rmsd)
          fastas.extend(fst)


f = open('afdb_ploop_final.csv', 'w')
data = list()
data.append(names)
data.append(queries)
data.append(seqs)
data.append(inits)
data.append(ends)
data.append(bfactors)
data.append(gkts)
data.append(meanbfactors)
data.append(IDs)
data.append(rmsds)
data.append(fastas)
data.append(jumps)
data.append(PorAs)
data.append(dires)
data.append(SheetIDs)
data.append(all_paras)
data.append(all_antis)
data.append(cycles)
data.append(undires)
data.append(withbranchs)
data.append(consecutives)
data.append(nstrandss)

data2 = list(zip(*data))
writer = csv.writer(f)
writer.writerows(data2)
f.close()
