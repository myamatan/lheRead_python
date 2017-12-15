import sys, math

import pylhe
from array import array

from ROOT import TH1D, gROOT, std
from ROOT import TFile, TTree, TCanvas
from ROOT import TLorentzVector, TVector3, TLorentzRotation

from sympy.functions.special.delta_functions import Heaviside

def printYogen(vec):
    print(vec.Px(), vec.Py(), vec.Pz(), vec.E())

inputDire='input/'
outputDire='output/'

if __name__ == "__main__":
    
    args = sys.argv
    f = TFile(outputDire+args[1]+'.root', 'recreate')
    tr = TTree('coll','coll')
    
    Mll = array("f", [0])
    Mjj = array("f", [0])
    Mlljj = array("f", [0])
    theta1 = array("f", [0])
    theta2 = array("f", [0])
    phi = array("f", [0])
    phi1 = array("f", [0])
    thetaStar = array("f", [0])
    MinPtRatio = array("f", [0])
    SqrtSSPt = array("f", [0])
    nbtag = array("f", [0])
    nbtagInFatjet = array("f", [0])
    ftag1 = array("f", [0])
    ftag2 = array("f", [0])

    tr.Branch('Mll', Mll, 'Mll/F')
    tr.Branch('Mjj', Mjj, 'Mjj/F')
    tr.Branch('theta1', theta1, 'theta1/F')
    tr.Branch('theta2', theta2, 'theta2/F')
    tr.Branch('phi', phi, 'phi/F')
    tr.Branch('phi1', phi1, 'phi1/F')
    tr.Branch('thetaStar', thetaStar, 'thetaStar/F')
    tr.Branch('Mlljj', Mlljj, 'Mlljj/F')
    tr.Branch('MinPtRatio', MinPtRatio, 'MinPtRatio/F')
    tr.Branch('SqrtSSPt', SqrtSSPt, 'SqrtSSPt/F')
    tr.Branch('nbtag', nbtag, 'nbtag/F')
    tr.Branch('nbtagInFatjet', nbtagInFatjet, 'nbtagInFatjet/F')
    tr.Branch('ftag1', ftag1, 'ftag1/F')
    tr.Branch('ftag2', ftag2, 'ftag2/F')

    print('Processing "',inputDire+args[1]+'.lhe', '....')

    if args[1].find('Zjets')!=-1 : isbkg=True
    else : isbkg=False

    for e in pylhe.readLHE (inputDire+args[1]+'.lhe'):
        nb = 0
        X = TLorentzVector(0.,0.,0.,0.)
        Zll = TLorentzVector(0.,0.,0.,0.)
        Zjj = TLorentzVector(0.,0.,0.,0.)
        l1 = TLorentzVector(0.,0.,0.,0.)
        l2 = TLorentzVector(0.,0.,0.,0.)
        q1 = TLorentzVector(0.,0.,0.,0.)
        q2 = TLorentzVector(0.,0.,0.,0.)

        isMuon = False
        bosonCount = 0
        for i in range( len(e['particles']) ):
            # Check if final state
            if e['particles'][i]['status'] == 1:
                pdgid = e['particles'][i]['id']
                Px =  e['particles'][i]['px']
                Py = e['particles'][i]['py']
                Pz = e['particles'][i]['pz']
                E = e['particles'][i]['e']
                if pdgid == -13 or pdgid == -11:
                    l1.SetPx(Px)
                    l1.SetPy(Py)
                    l1.SetPz(Pz)
                    l1.SetE(E)
                    if pdgid == -11:
                        isMuon = False
                    else:
                        isMuon = True
                elif pdgid == 13 or pdgid == 11:
                    l2.SetPx(Px)
                    l2.SetPy(Py)
                    l2.SetPz(Pz)
                    l2.SetE(E)
                elif isbkg==False:
                    if abs(pdgid) < 7 and pdgid < 0:
                        q1.SetPx(Px)
                        q1.SetPy(Py)
                        q1.SetPz(Pz)
                        q1.SetE(E)
                        ft1 = abs(pdgid)
                    elif abs(pdgid) < 7 and pdgid > 0:
                        q2.SetPx(Px)
                        q2.SetPy(Py)
                        q2.SetPz(Pz)
                        q2.SetE(E)
                        ft2 = abs(pdgid)
                else :
                    if pdgid == 21 or (abs(pdgid) < 7 and pdgid < 0):
                        q1.SetPx(Px)
                        q1.SetPy(Py)
                        q1.SetPz(Pz)
                        q1.SetE(E)
                        ft1 = abs(pdgid)
                    elif pdgid == 21 or (abs(pdgid) < 7 and pdgid > 0):
                        q2.SetPx(Px)
                        q2.SetPy(Py)
                        q2.SetPz(Pz)
                        q2.SetE(E)
                        ft2 = abs(pdgid)

                if abs(pdgid) == 5:
                    nb+=1
        
        j1= TLorentzVector(0.,0.,0.,0.)
        j2= TLorentzVector(0.,0.,0.,0.)
        j1 = q1
        j2 = q2
        if q1.Pt() > q2.Pt():
            q1 = j1
            q2 = j2
            ftag1[0] = ft1
            ftag2[0] = ft2
        else:
            q1 = j2
            q2 = j1
            ftag1[0] = ft2
            ftag2[0] = ft1

        Zll = l1 + l2
        Zjj = q1 + q2
        X = Zll + Zjj
        
        # Gev to MeV
        X *= 1e+3
        Zll *= 1e+3
        Zjj *= 1e+3
        l1 *= 1e+3
        l2 *= 1e+3
        q1 *= 1e+3
        q2 *= 1e+3

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
       
        if q1.DeltaR(q2) > 1.0:
            G *= Heaviside(4.5-abs(q1.Eta()))
            G *= Heaviside(4.5-abs(q2.Eta()))
            G *= Heaviside(q1.Pt()-30.0e+3)
            G *= Heaviside(q2.Pt()-30.0e+3)
        else:
            G *= Heaviside( (q1+q2).Pt()-200.0e+3 )
            G *= Heaviside( 2.0-abs((q1+q2).Eta()) )

        
        # Event selection
        Es = 1.0
        if isMuon:
            Es *= Heaviside(99e+3-(l1+l2).M())
            Es *= Heaviside((l1+l2).M()-83e+3)
        else:
            Es *= Heaviside(94e+3+0.01850*(l1+l2).Pt()-(l1+l2).M())
            Es *= Heaviside((l1+l2).M()-85.63e+3+0.01170*(l1+l2).Pt())
        if q1.DeltaR(q2) > 1.0:
            if q1.Pt() > q2.Pt():
                Es *= Heaviside(q1.Pt()-60e+3)
            else:
                Es *= Heaviside(q2.Pt()-60e+3)

        if q1.DeltaR(q2) > 1.0:
            MinPtRatio[0] = -99.
            SqrtSSPt[0] = math.sqrt((l1+l2).Pt()**2+(q1+q2).Pt()**2)/(l1+l2+q1+q2).M()
            nbtag[0] = nb
            nbtagInFatjet[0] = -99
        else:
            MinPtRatio[0] = min( (l1+l2).Pt(), (q1+q2).Pt() ) / (l1+l2+q1+q2).M() 
            SqrtSSPt[0] = -99.
            nbtag[0] = nb
            nbtagInFatjet[0] = nb

        #Es *= Heaviside(np.sqrt((l1+l2).Pt()**2+(q1+q2).Pt()**2)/(l1+l2+q1+q2).M()-0.4)
        #Es *= Heaviside(105e+3-(q1+q2).M())
        #Es *= Heaviside((q1+q2).M()-70e+3)
        
        # Invariant mass of Z1 and Z2
        Mll[0] = (l1 + l2).M()
        Mjj[0] = (q1 + q2).M()
        Mlljj[0] = (l1 + l2 + q1 + q2).M()

        # Boost to X's static system
        l_X = TLorentzRotation().Boost(-1.*X.Px()/X.E(), -1.*X.Py()/X.E(), -1.*X.Pz()/X.E())
        X = l_X * X;
        Z1 = l_X * Zll;
        Z2 = l_X * Zjj;
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
       
        # Five angular information (theta1, theta2, Phi, Phi1, thetaStar)
        theta1[0] = math.acos((-1.)*Z2_Z1SS.Vect().Unit().Dot(l1_Z1SS.Vect().Unit()))

        theta2[0] = math.acos((-1.)*Z1_Z2SS.Vect().Unit().Dot(q1_Z2SS.Vect().Unit()))
        
        v1ortho = l1.Vect().Cross(l2.Vect())
        v2ortho = q1.Vect().Cross(q2.Vect())
        n1 = v1ortho.Unit()
        n2 = v2ortho.Unit()
        A = Z1.Vect().Dot(n1.Cross(n2))
        if A!=0 : phi[0] = A/abs(A)*math.acos((-1.)*n1.Dot(n2))
        else : phi[0] = -99
        
        nz = TVector3(0,0,1)
        nsc = (nz.Cross(Z1.Vect())).Unit()
        B = Z1.Vect().Dot(n1.Cross(nsc))
        phi1[0] = B/abs(B)*math.acos(n1.Dot(nsc))
        
        thetaStar[0] = Z1.Vect().Theta()


        #if G > 0 and Es > 0:
        if G > 0:
            tr.Fill()
            
    f.Write()
    f.Close()
