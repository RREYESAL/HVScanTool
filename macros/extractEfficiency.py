#!/usr/bin/env python
from ROOT import *
import re
import os, time

runsToAnalyze = [272818]

## Build bidirectional mapping between detId and roll name
detIdToRollName, rollNameToDetId = {}, {}
for l in open("data/detId.txt").readlines():
    l = l.strip()
    if len(l) == 0 or l[0] == '#': continue
    l = l.split()
    if len(l) != 2: continue

    rollName, detId = l[0], int(l[1])
    detIdToRollName[detId] = rollName
    rollNameToDetId[rollName] = detId

## Read HV setups
hvSetup = {}
for l in open("data/hvEffective.txt").readlines():
    l = l.strip()
    if len(l) == 0 or l[0] == '#': continue
    l = l.split()
    if len(l) != 8: continue

    fill, run, lsBegin, lsEnd = [int(x) for x in l[:4]]
    hvEff, hvSet, pressure = [float(x) for x in l[5:]]
    rollNamePattern = l[4].replace('.','?').replace('*', '.*')

    if run not in runsToAnalyze: continue
    hvSetup[(fill, run, lsBegin, lsEnd)] = [rollNamePattern, hvEff, hvSet, pressure]

data = {}

## Analyze data
for l in open("data/inputFiles.txt").readlines():
    l = l.strip()
    if len(l) == 0 or l[0] == '#': continue
    l = l.split()
    if len(l) != 5: continue

    fill, run, lsBegin, lsEnd = [int(x) for x in l[:4]]
    fileName = l[4]

    if run not in runsToAnalyze: continue
    if not os.path.exists('/eos/cms/%s' % fileName): continue

    if (fill, run, lsBegin, lsEnd) not in hvSetup: continue
    rollNamePattern, hvEff, hvSet, pressure = hvSetup[(fill, run, lsBegin, lsEnd)]

    f = TFile('/eos/cms/%s' % fileName)
    tB = f.Get("effTreeBarrel")
    tE = f.Get("effTreeEndcap")

    for entry in tB:
        eff, effErr = entry.fiducialCutEff, entry.fiducialCutEffErr
        detId = entry.detId
        if detId not in detIdToRollName:
            print detId, ' not in the DB'
            continue
        detName = detIdToRollName[detId]
        if not re.match(rollNamePattern, detName): continue
        if detName not in data: data[detName] = []
        data[detName].append([eff, effErr, hvEff])

    for entry in tE:
        eff, effErr = entry.fiducialCutEff, entry.fiducialCutEffErr
        detId = entry.detId
        if detId not in detIdToRollName:
            print detId, ' not in the DB'
            continue
        detName = detIdToRollName[detId]
        if not re.match(rollNamePattern, detName): continue
        if detName not in data: data[detName] = []
        data[detName].append([eff, effErr, hvEff])

if not os.path.exists("res"): os.mkdir("res")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
c = TCanvas("c", "c", 600, 400)
hFrame = TH1F("hFrame", "hFrame", 100, 8500, 10000)
hFrame.SetMinimum(0)
hFrame.SetMaximum(110)
grps = {}
sigmoid = "[0]/(1.0+TMath::Exp([1]*(x-[2])))"
ftn = TF1("ftn", sigmoid, 8500, 10000)
ftn.SetParNames("emax", "slope", "hv50")
for d in data:
    g = TGraphErrors()
    for e, err, v in data[d]:
        n = g.GetN()
        g.SetPoint(n, v, 100*e)
        g.SetPointError(n, 0, 100*err)
    ftn.SetParLimits(0, 0.0, 99.999) #bound emax parameter 
    ftn.SetParameters(95, -100./(10000-8500), 9000)

    g.Fit(ftn, "BR")
    g.Fit(ftn, "BR")
    res = g.Fit(ftn, "BRMS")

    hFrame.Draw()
    g.Draw("same,P")

    c.Print("res/fit_%s.png" % d)

    grps[d] = [g, res]

    time.sleep(0.5)

