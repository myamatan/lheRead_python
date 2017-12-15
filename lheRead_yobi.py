import sys, math

import pylhe
import numpy as np
from array import array

from ROOT import TH1D, gROOT, std
from ROOT import TFile, TTree, TCanvas
from ROOT import TLorentzVector, TVector3, TLorentzRotation

from sympy.functions.special.delta_functions import Heaviside

def printYogen(vec):
    print(vec.Px(), vec.Py(), vec.Pz(), vec.E())

inputDire='input_mg5/'
outputDire='output_mg5/'

if __name__ == "__main__":
    
    args = sys.argv
    f = TFile(outputDire+args[1]+'.root', 'recreate')
    tr = TTree('coll','coll')
    
    mll = array("f", [0])
    mqq = array("f", [0])
    theta1 = array("f", [0])
    theta2 = array("f", [0])
    Phi = array("f", [0])
    Phi1 = array("f", [0])
    thetaStar = array("f", [0])

    tr.Branch('mll', mll, 'mll/F')
    tr.Branch('mqq', mqq, 'mqq/F')
    tr.Branch('theta1', theta1, 'theta1/F')
    tr.Branch('theta2', theta2, 'theta2/F')
    tr.Branch('Phi', Phi, 'Phi/F')
    tr.Branch('Phi1', Phi1, 'Phi1/F')
    tr.Branch('thetaStar', thetaStar, 'thetaStar/F')

    print('Processing "',inputDire+args[1]+'.lhe', '....')

    for e in pylhe.readLHE (inputDire+args[1]+'.lhe'):
        X = TLorentzVector(0.,0.,0.,0.)
        Z1 = TLorentzVector(0.,0.,0.,0.)
        Z2= TLorentzVector(0.,0.,0.,0.)
        l1= TLorentzVector(0.,0.,0.,0.)
        l2= TLorentzVector(0.,0.,0.,0.)
        q1= TLorentzVector(0.,0.,0.,0.)
        q2= TLorentzVector(0.,0.,0.,0.)
        ID1 = 0
        ID2 = 0
        isMuon = False
        bosonCount = 0
        for i in range( len(e['particles']) ):
            if e['particles'][i]['id'] == -13 or e['particles'][i]['id'] == -11:
                l1.SetPx(e['particles'][i]['px'])
                l1.SetPy(e['particles'][i]['py'])
                l1.SetPz(e['particles'][i]['pz'])
                l1.SetE(e['particles'][i]['e'])
                ID1 = int(e['particles'][i]['mother1'] - 1)
                if e['particles'][i]['id'] == -11:
                    isMuon = False
                else:
                    isMuon = True 
            elif e['particles'][i]['id'] == 13 or e['particles'][i]['id'] == 11:
                l2.SetPx(e['particles'][i]['px'])
                l2.SetPy(e['particles'][i]['py'])
                l2.SetPz(e['particles'][i]['pz'])
                l2.SetE(e['particles'][i]['e'])
            elif abs(e['particles'][i]['id']) < 7 and e['particles'][i]['id'] < 0:
                q1.SetPx(e['particles'][i]['px'])
                q1.SetPy(e['particles'][i]['py'])
                q1.SetPz(e['particles'][i]['pz'])
                q1.SetE(e['particles'][i]['e'])
                ID2 = int(e['particles'][i]['mother1'] - 1)
            elif abs(e['particles'][i]['id']) < 7 and e['particles'][i]['id'] > 0:
                q2.SetPx(e['particles'][i]['px'])
                q2.SetPy(e['particles'][i]['py'])
                q2.SetPz(e['particles'][i]['pz'])
                q2.SetE(e['particles'][i]['e'])    

        Z1.SetPx(e['particles'][ID1]['px'])
        Z1.SetPy(e['particles'][ID1]['py'])
        Z1.SetPz(e['particles'][ID1]['pz'])
        Z1.SetE(e['particles'][ID1]['e'])
        Z2.SetPx(e['particles'][ID2]['px'])
        Z2.SetPy(e['particles'][ID2]['py'])
        Z2.SetPz(e['particles'][ID2]['pz'])
        Z2.SetE(e['particles'][ID2]['e'])
        
        X = Z1 + Z2

        # Acceptance cut
        #G = Heaviside(2.5-abs(l1.Eta()))
        #G *= Heaviside(2.5-abs(l2.Eta()))
        #G *= Heaviside(2.5-abs(q1.Eta()))
        #G *= Heaviside(2.5-abs(q2.Eta()))
        
        # Gev to MeV
        X *= 1e+3
        Z1 *= 1e+3
        Z2 *= 1e+3
        l1 *= 1e+3
        l2 *= 1e+3
        q1 *= 1e+3
        q2 *= 1e+3
        
        # Boost to X's static system
        l_X = TLorentzRotation().Boost(-1.*X.Px()/X.E(), -1.*X.Py()/X.E(), -1.*X.Pz()/X.E())
        X = l_X * X;
        Z1 = l_X * Z1;
        Z2 = l_X * Z2;
        l1 = l_X * l1;
        l2 = l_X * l2;
        q1 = l_X * q1;
        q2 = l_X * q2;
    
        # Boost to Z1's static system
        l_Z1 = TLorentzRotation().Boost(-1.*Z1.Px()/Z1.E(), -1.*Z1.Py()/Z1.E(), -1.*Z1.Pz()/Z1.E())
        Z1_Z1SS = l_Z1 * Z1
        Z2_Z1SS = l_Z1 * Z2
        l1_Z1SS = l_Z1 * l1
        l2_Z1SS = l_Z1 * l2
        
        # Boost to Z2's static system
        l_Z2 = TLorentzRotation().Boost(-1.*Z2.Px()/Z2.E(), -1.*Z2.Py()/Z2.E(), -1.*Z2.Pz()/Z2.E())
        Z1_Z2SS = l_Z2 * Z1
        Z2_Z2SS = l_Z2 * Z2
        q1_Z2SS = l_Z2 * q1
        q2_Z2SS = l_Z2 * q2
       
        # Invariant mass of Z1 and Z2
        mll[0] = (l1 + l2).M()
        mqq[0] = (q1 + q2).M()
        
        # Five angular information (theta1, theta2, Phi, Phi1, thetaStar)
        theta1[0] = math.acos((-1.)*Z2_Z1SS.Vect().Unit().Dot(l1_Z1SS.Vect().Unit()))
        theta1p = math.acos((-1.)*Z2_Z1SS.Vect().Unit().Dot(l2_Z1SS.Vect().Unit()))

        theta2[0] = math.acos((-1.)*Z1_Z2SS.Vect().Unit().Dot(q1_Z2SS.Vect().Unit()))
        theta2p = math.acos((-1.)*Z1_Z2SS.Vect().Unit().Dot(q2_Z2SS.Vect().Unit()))
        
        v1ortho = l1.Vect().Cross(l2.Vect())
        v2ortho = q1.Vect().Cross(q2.Vect())
        n1 = v1ortho.Unit()
        n2 = v2ortho.Unit()
        A = Z1.Vect().Dot(n1.Cross(n2))
        Phi[0] = A/abs(A)*math.acos((-1.)*n1.Dot(n2))
        
        nz = TVector3(0,0,1)
        nsc = (nz.Cross(Z1.Vect())).Unit()
        B = Z1.Vect().Dot(n1.Cross(nsc))
        Phi1[0] = B/abs(B)*math.acos(n1.Dot(nsc))
        
        thetaStar[0] = Z1.Vect().Theta()
        thetaStar2 = Z2.Vect().Theta()

        # Objects definition cuts
        G = 1.0
        if isMuon:
            G *= Heaviside(2.5-abs(l1.Eta()))
            G *= Heaviside(2.5-abs(l2.Eta()))
            G *= Heaviside(l1.Pt()-7.0e+3)
            G *= Heaviside(l2.Pt()-7.0e+3)
        else:
            G *= Heaviside(2.47-abs(l1.Eta()))
            G *= Heaviside(2.47-abs(l2.Eta()))
            G *= Heaviside(l1.Pt()-7.0e+3)
            G *= Heaviside(l2.Pt()-7.0e+3)
        
        G *= Heaviside(4.5-abs(q1.Eta()))
        G *= Heaviside(4.5-abs(q2.Eta()))
        G *= Heaviside(q1.Pt()-30.0e+3)
        G *= Heaviside(q2.Pt()-30.0e+3)
        
        if G > 0:
            tr.Fill()

    f.Write()
    f.Close()
        
