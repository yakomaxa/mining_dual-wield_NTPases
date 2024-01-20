import foldcomp
import pymol
from pymol import stored
import numpy as np
import re
import csv
import sys
args = sys.argv

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
                IDs.append(name)
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


#/home/ksakuma/myprojects/FOLDCOMP/database/afdb_foldcomp

#idlist=open(args[1],mode="r")
#ids=idlist.read().splitlines()
#idlist.close()

with foldcomp.open("./splitdb/afdb_uniprot_v4") as db:
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
  count = 0
  pcount =0
  for (name_tmp, pdb) in db:
      name=name_tmp.replace(",","").replace("/","").replace(" ","_")
      abego,fasta,b=pdb2abego(pdb)
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

      for qabego in ["BBBBBBGAGAAAAA","BBBEBBGAGAAAAA","BBBABBGAGAAAAA","BBBGBBGAGAAAAA"]:
          myname, query, seq, init, end, bfactor, gkt, meanbfactor, ID, rmsd, fst = iterative_search(name=name,fasta=fasta,abego=abego,query=qabego,bfactor=b)

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
data2 = list(zip(*data))
writer = csv.writer(f)
writer.writerows(data2)
f.close()
