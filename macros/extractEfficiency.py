#!/usr/bin/env python
from ROOT import *
import re
import os, time

#runsToAnalyze = [317591,317626] ## HVScan A 2018
runsToAnalyze = [303948] ## HVScan B 2017
#runsToAnalyze = [295602, 295603, 295605] ## HVScan A 2017
#runsToAnalyze = [272818] ## HVScan 2016

## Build bidirectional mapping between detId and roll name
detIdToDetName, detNameToDetId = {}, {}
for l in open("data/detId.txt").readlines():
    l = l.rstrip('#').strip()
    if len(l) == 0 or l[0] == '#': continue
    l = l.split()
    if len(l) != 2: continue

    detName, detId = l[0], int(l[1])
    detIdToDetName[detId] = detName
    detNameToDetId[detName] = detId

## Read HV setups
hvSetup = {}
for l in open("data/hvEffective.txt").readlines():
    l = l.rstrip('#').strip()
    if len(l) == 0 or l[0] == '#': continue
    l = l.split()
    if len(l) != 8: continue

    fill, run, lsBegin, lsEnd = [int(x) for x in l[:4]]
    hvEff, hvSet, pressure = [float(x) for x in l[5:]]
    detNamePattern = l[4].replace('.','?').replace('*', '.*')

    if run not in runsToAnalyze: continue
    hvSetup[(fill, run, lsBegin, lsEnd)] = [detNamePattern, hvEff, hvSet, pressure]

data = {}

## Analyze data
for l in open("data/inputFiles.txt").readlines():
    l = l.rstrip('#').strip()
    if len(l) == 0 or l[0] == '#': continue
    l = l.split()
    if len(l) != 5: continue

    fill, run, lsBegin, lsEnd = [int(x) for x in l[:4]]
    fileName = l[4]

    if run not in runsToAnalyze: continue
    if not os.path.exists('/eos/cms/%s' % fileName): continue

    if (fill, run, lsBegin, lsEnd) not in hvSetup: continue
    detNamePattern, hvEff, hvSet, pressure = hvSetup[(fill, run, lsBegin, lsEnd)]

    f = TFile('/eos/cms/%s' % fileName)
    tB = f.Get("effTreeBarrel")
    tE = f.Get("effTreeEndcap")

    for entry in tB:
        eff, effErr = entry.fiducialCutEff, entry.fiducialCutEffErr
        cls = entry.clustersize
        if eff == 0 or cls == 0: continue
        detId = entry.detId
        if detId not in detIdToDetName:
            print detId, ' not in the DB'
            continue
        detName = detIdToDetName[detId]
        if not re.match(detNamePattern, detName): continue
        if detName not in data: data[detName] = []
        data[detName].append([eff, effErr, cls, hvEff])

    for entry in tE:
        eff, effErr = entry.fiducialCutEff, entry.fiducialCutEffErr
        cls = entry.clustersize
        if eff == 0 or cls == 0: continue
        detId = entry.detId
        if detId not in detIdToDetName:
            print detId, ' not in the DB'
            continue
        detName = detIdToDetName[detId]
        if not re.match(detNamePattern, detName): continue
        if detName not in data: data[detName] = []
        data[detName].append([eff, effErr, cls, hvEff])

if not os.path.exists("res"): os.mkdir("res")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

html = open("res/index.html", "w")
print>>html, """<html>
<head><title>HVScan analysis fit results</title></head>
<body>
"""

c = TCanvas("c", "c", 600, 400)
hFrame = TH1F("hFrame", "hFrame;HV [V];Efficiency [%]", 100, 8400, 10000)
hFrame.SetMinimum(0)
hFrame.SetMaximum(110)

cBads = {}
fitResults = {}

hvMin, hvMax, clsMax, effMax = 8400, 10000, 10, 110

sigmoid = "[0]/(1.0+TMath::Exp([1]*(x-[2])))"
ftn = TF1("ftn", sigmoid, hvMin, hvMax)
ftn.SetParNames("emax", "slope", "hv50")

clsScale = float(clsMax)/effMax
clsAxis = TGaxis(hvMax, 0, hvMax, effMax, 0, clsMax, 510,"+L");
clsAxis.SetLineColor(kBlue)
clsAxis.SetLabelColor(kBlue)
clsAxis.SetTitle("Cluster size")

for detName in data:
    if len(data[detName]) < 4: continue

    gEff = TGraphErrors()
    gCls = TGraph()
    for e, err, cls, v in data[detName]:
        n = gEff.GetN()
        gEff.SetPoint(n, v, 100*e)
        gEff.SetPointError(n, 0, 100*err)
        gCls.SetPoint(n, v, cls/clsScale)
    ftn.SetParLimits(0, 0.0, 99.999) #bound emax parameter 
    ftn.SetParLimits(1, -0.02, 0) #bound slope parameter 
    ftn.SetParLimits(2, hvMin, hvMax)
    ftn.SetParameters(99, -1e-9, 9200)

    gEff.Fit(ftn, "QBRNEX0")
    gEff.Fit(ftn, "QBRNEX0")
    fitResult = gEff.Fit(ftn, "QBRMSNEX0")

    fitResults[detName] = [fitResult]

    #time.sleep(0.5)

    gCls.SetLineColor(kBlue)
    gCls.SetMarkerColor(kBlue)
    gCls.SetMarkerStyle(kFullSquare)
    gCls.SetMarkerSize(0.5)

    c.cd()
    hFrame.Draw()
    gEff.Draw("same,P")
    gCls.Draw("same,LP")
    clsAxis.Draw()
    ftn.Draw("same")

    tl = TPaveText(9400,60,9800, 80, "")
    tl.SetBorderSize(0)
    tl.SetFillStyle(0)
    tl.AddText(detName)
    tl.AddText("#epsilon_{max} = %.0f#pm%.1f" % (fitResult.Value(0), fitResult.Error(0)))
    tl.AddText("slope = %.3f#pm%.4f" % (fitResult.Value(1), fitResult.Error(1)))
    tl.AddText("HV_{50} = %.0f#pm%.1f" % (fitResult.Value(2), fitResult.Error(2)))
    tl.Draw()

    c.Print("res/fit_%s.png" % detName)

    if not fitResult.IsValid():
        cBad = TCanvas("cBad_%s" % detName, "Bad fit result %s" % detName, 600, 400)
        hFrame.Draw()
        gBad = gEff.Clone()
        gClsBad = gBad.Clone()
        ftnBad = ftn.Clone("%s_%s" % (ftn.GetName(), detName))
        gBad.Draw("same,P")
        gClsBad.Draw("same,LP")
        ftnBad.Draw("same")
        clsAxis.Draw()
        tl.Draw()
        cBads[detName] = [cBad, hFrame, gBad, gClsBad, ftnBad, tl]

    if fitResult.IsValid():
        print>>html, ('<a href="fit_%s.png"><img src="fit_%s.png" style="width:300px" /></a>' % (detName, detName))
    else:
        print>>html, ('<a href="fit_%s.png"><img src="fit_%s.png" style="width:300px;border:1px solid red" /></a>' % (detName, detName))

    #break

hEmaxRB = TH1D("hEmaxRB", "Barrel E_{max};Efficiency at plataeu [%%];Number of rolls", 51, 50, 101)
hEmaxRE = TH1D("hEmaxRE", "Endcap E_{max};Efficiency at plataeu [%%];Number of rolls", 51, 50, 101)
hEmaxRE4 = TH1D("hEmaxRE4", "Endcap disk4 E_{max};Efficiency at plataeu [%%];Number of rolls", 51, 50, 101)
hEmaxRB.SetFillColor(kRed+1)
hEmaxRE.SetFillColor(kBlue+2)
hEmaxRE4.SetFillColor(kGreen+3)
for detName in fitResults:
    if detName in cBads: continue
    fitResult = fitResults[detName]
    if len(fitResult) == 0: continue # this should never happend
    fitResult = fitResult[0]

    emax, emaxErr = fitResult.Value(0), fitResult.Error(0)
    slope, slopeErr = fitResult.Value(1), fitResult.Error(1)
    hv50, hv50Err = fitResult.Value(2), fitResult.Error(2)

    if detName.startswith("W"): hEmaxRB.Fill(emax)
    elif detName.startswith("RE+4") or detName.startswith("RE-4"): hEmaxRE4.Fill(emax)
    else: hEmaxRE.Fill(emax)
hsEmax = THStack("hsEmax", "hsEmax")
hsEmax.Add(hEmaxRB)
hsEmax.Add(hEmaxRE)
hsEmax.Add(hEmaxRE4)
cEmax = TCanvas("cEmax", "cEmax", 500, 500)
hsEmax.Draw()

print>>html, """</body>
</html>
"""
