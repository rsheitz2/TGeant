#!/usr/bin/python
import sys
import ROOT

'''
Used to get cut statistics
root files locations should be listed in a .txt file

Warning:
At the moment, must go back root 5.34.24/x86_64-slc6-gcc47-opt and 
and compiler 4.7.2/x86_64-slc6-gcc47-opt
execute ipythonSetup.sh to do this
Usage:
ipython 
%run CutStats_FilesFromTxt.py PhastData.txt UserEvent# MainData.txt
'''
if len(sys.argv) < 3:
    print "Used to get cut statistics"
    print "root files locations should be listed in a .txt file in order from",
    print "W07 to W15"
    print " "
    print "Warning:"
    print "At the moment, must go back root 5.34.24/x86_64-slc6-gcc47-opt and"
    print "and compiler 4.7.2/x86_64-slc6-gcc47-opt"
    print "execute ipythonSetup.sh to do this"
    print " "
    print "Enter text file with names of Phast root files to look through"
    print "Enter UserEvent number (i.e. 420)"
    print "Enter text file with names of Main files to look through"
else:
    #Open files
    phastFiles = open(sys.argv[1] )
    mainFiles = open(sys.argv[3] )

    #Setup labels and periods
    WLabels = []
    W07 = []
    W08 = []
    W09 = []
    W10 = []
    W11 = []
    W12 = []
    W13 = []
    W14 = []
    W15 = []
    WAll = []
    AllPeriods = [WLabels, W07, W08, W09, W10, W11, W12, W13, W14, W15, WAll]

    #Phast Files
    start_bin = 1
    bin_space = 10
    first = True
    file_iter = 1
    for line in phastFiles:
        separate = line.split("\n")
        oneFile = ROOT.TFile(separate[0] )
        user = "UserEvent" + sys.argv[2]
        
        hPhastCuts = oneFile.Get(user+"/DiMuonCuts")
        for bin in range(start_bin, hPhastCuts.GetXaxis().GetNbins()+1, bin_space):
            if (first):
                AllPeriods[0].append(hPhastCuts.GetXaxis().GetBinLabel(bin) )
            
            AllPeriods[file_iter].append(hPhastCuts.GetBinContent(bin) )

        first = False
        file_iter += 1

    #Main Files
    start_bin = 1
    bin_space = 10
    nCuts = 6
    first = True
    file_iter = 1
    for line in mainFiles:
        separate = line.split("\n")
        oneFile = ROOT.TFile(separate[0] )
                
        hCuts = oneFile.Get("hCuts")
        for bin in range(start_bin, nCuts*bin_space, bin_space):
            if (first):
                AllPeriods[0].append(hCuts.GetXaxis().GetBinLabel(bin+1) )
            AllPeriods[file_iter].append(hCuts.GetBinContent(bin) )

        first = False
        file_iter += 1

    #Print output to file
    fPeriodCuts = open('PeriodCuts.txt', 'w')
    for period in AllPeriods:
        for cut in period:
            fPeriodCuts.write(str(cut) + " ")

        fPeriodCuts.write("\n")
