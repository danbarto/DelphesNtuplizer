#!/usr/bin/env python
import sys
import ROOT
from collections import OrderedDict
from ROOT import TLorentzVector
from array import array
import numpy as np
import argparse


class TreeProducer:
    def __init__(self, debug):

         # flat tree branches
         self.debug = debug

         self.t = ROOT.TTree( "mytree","TestTree" )
         self.maxn = 9999

         self.fatjet_size                   = array('i', [0 ])
         self.fatjet_pt                     = array('f', self.maxn*[0. ])
         self.fatjet_eta                    = array('f', self.maxn*[0. ])
         self.fatjet_phi                    = array('f', self.maxn*[0. ])
         self.fatjet_mass                   = array('f', self.maxn*[0. ])
         self.fatjet_tau1                   = array('f', self.maxn*[0. ])
         self.fatjet_tau2                   = array('f', self.maxn*[0. ])
         self.fatjet_tau3                   = array('f', self.maxn*[0. ])
         self.fatjet_tau4                   = array('f', self.maxn*[0. ])
         self.fatjet_msoftdrop              = array('f', self.maxn*[0. ])

         self.pfcand_size                   = array('i', [0 ])
         self.pfcand_pt                     = array('f', self.maxn*[0. ])
         self.pfcand_eta                    = array('f', self.maxn*[0. ])
         self.pfcand_phi                    = array('f', self.maxn*[0. ])

         self.t.Branch("fatjet_size",   self.fatjet_size,   "fatjet_size/I")
         self.t.Branch("fatjet_pt", self.fatjet_pt, "fatjet_pt[fatjet_size]/F")
         self.t.Branch("fatjet_eta", self.fatjet_eta, "fatjet_eta[fatjet_size]/F")
         self.t.Branch("fatjet_phi", self.fatjet_phi, "fatjet_phi[fatjet_size]/F")
         self.t.Branch("fatjet_mass", self.fatjet_mass, "fatjet_mass[fatjet_size]/F")
         self.t.Branch("fatjet_pt", self.fatjet_pt, "fatjet_pt[fatjet_size]/F")
         self.t.Branch("fatjet_tau1", self.fatjet_tau1, "fatjet_tau1[fatjet_size]/F")
         self.t.Branch("fatjet_tau2", self.fatjet_tau2, "fatjet_tau2[fatjet_size]/F")
         self.t.Branch("fatjet_tau3", self.fatjet_tau3, "fatjet_tau3[fatjet_size]/F")
         self.t.Branch("fatjet_tau4", self.fatjet_tau4, "fatjet_tau4[fatjet_size]/F")
         self.t.Branch("fatjet_msoftdrop", self.fatjet_msoftdrop, "fatjet_msoftdrop[fatjet_size]/F")

         self.t.Branch("pfcand_size",   self.pfcand_size,   "pfcand_size/I")
         self.t.Branch("pfcand_pt",     self.pfcand_pt,     "pfcand_pt[pfcand_size]/F")
         self.t.Branch("pfcand_eta",    self.pfcand_eta,    "pfcand_eta[pfcand_size]/F")
         self.t.Branch("pfcand_phi",    self.pfcand_phi,    "pfcand_phi[pfcand_size]/F")


    def fill(self):
        self.t.Fill()

    def write(self):
        self.t.Write()

#_______________________________________________________
def dr_match(p1, p2, drmin):
    dr = p1.P4().DeltaR(p2.P4())
    return dr < drmin


#_____________________________________________________________________________________________________________
def main():

    ROOT.gSystem.Load("libDelphes")
    try:
      ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
      ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
    except:
      pass

    parser = argparse.ArgumentParser()
    parser.add_argument ('-i', '--input', help='input Delphes files',  default='delphes.root')
    parser.add_argument ('-o', '--output', help='output flat tree',  default='tree.root')
    parser.add_argument ('-n', '--nev', help='number of events', type=int, default=-1)
    parser.add_argument ('-d', '--debug', help='debug flag',  action='store_true',  default=False)

    args = parser.parse_args()

    inputFile = args.input
    outputFile = args.output
    nevents = args.nev
    debug = args.debug

    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)
    

    # Create object of class ExRootTreeReader
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    ## for now only M for electrons, LT for muons and LT for photons are defined !!
    ## should dervie new parameterisations for other working points

    branchVertex          = treeReader.UseBranch('Vertex')   
    branchParticle        = treeReader.UseBranch('Particle') 
    branchGenJet          = treeReader.UseBranch('GenJet')   

    branchFatJet          = treeReader.UseBranch('JetPUPPIAK8') 

    # NEED these branches to access jet constituents
    branchPuppiCandidate  = treeReader.UseBranch('ParticleFlowCandidate')

    branchRho             = treeReader.UseBranch('Rho')

    treeProducer = TreeProducer(debug)

    if nevents > 0:
        numberOfEntries = nevents

    ################ Start event loop #######################
    for entry in range(0, numberOfEntries):

        # Load selected branches with data from specified event
        treeReader.ReadEntry(entry)

        if (entry+1)%1000 == 0:
            print ' ... processed {} events ...'.format(entry+1)

        for item in branchFatJet:

            i = 0   
            jetp4 = item.P4()
            softDrop = TLorentzVector() 
            softDrop = item.SoftDroppedJet
            treeProducer.fatjet_pt                     [i] = jetp4.Pt()
            treeProducer.fatjet_eta                    [i] = jetp4.Eta()
            treeProducer.fatjet_phi                    [i] = jetp4.Phi()
            treeProducer.fatjet_mass                   [i] = jetp4.M()
            treeProducer.fatjet_tau1                   [i] = item.Tau[0]            
            treeProducer.fatjet_tau2                   [i] = item.Tau[1]            
            treeProducer.fatjet_tau3                   [i] = item.Tau[2]            
            treeProducer.fatjet_tau4                   [i] = item.Tau[3]            
            treeProducer.fatjet_msoftdrop              [i] = softDrop.M()      

            treeProducer.fatjet_size[0] = 1


            for j in xrange(len(item.Constituents)):
                const = item.Constituents.At(j)
                p4 = ROOT.TLorentzVector(0., 0., 0., 0.)
                if isinstance(const, ROOT.ParticleFlowCandidate):
                    p4 = ROOT.ParticleFlowCandidate(const).P4()
                    #nconst +=1
                    treeProducer.pfcand_pt      [j] = p4.Pt()
                    treeProducer.pfcand_phi     [j] = p4.Phi()
                    treeProducer.pfcand_eta     [j] = p4.Eta()
                    if debug: print '       PFCandidate: ',const.PID, p4.Pt(), p4.Eta(), p4.Phi(), p4.M()

            treeProducer.pfcand_size[0] = j+1

            ## fill tree 
            treeProducer.fill()


    out_root = ROOT.TFile(outputFile,"RECREATE")
    out_root.cd()
    treeProducer.write()
 

#_______________________________________________________________________________________
if __name__ == "__main__":
    main()

