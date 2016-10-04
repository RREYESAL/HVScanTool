#HV Scan Tool

The HV Scan Tool it is a series of codes that you can run depending on the output you need. 

1. To start the analysis you need to get the input root files from: 

eos ls /store/group/dpg_rpc/comm_rpc/Run-II/data2016/HVScan_May9/ForHVscanAna

the name of files are: AnalyzeEfficiency_272818_RPCMon_p\*.root 
They are 12 in total, also you have to create a file a hvEffective.txt with the values of the Effective High Voltages in the taking data. The values must to be in kV  
You can get more information in the twiki page: https://twiki.cern.ch/twiki/bin/view/CMS/RPCHvscan2016
Put the input files in the data directory. 

2. Once you have the input files in the data directory go to macro directory, the next step is run the FitData.C macro. This macro needs a input parameter depending on where you want 
   to do the fit "barrel" or "endcap", for instance you can do in the prompt line $ root -l -q 'FitData.C("barrel")'. The output is the analysis of the data, this creates the results directory 
   and 1 subdirectory per chamber in the RPC with the information of the fit results of the efficiency and the cluster size in two txt Files "fitData.txt" and "fitDataCls.txt". 
   Basically the macro takes the root files, next it does two distributions, efficiency vs HV and clustersize vs HV, and then fit them using a Sigmoid function and a polynomial function 
   respectively. 

3. You can see the plots of the distributions and the fit running the pngProducer.C. This macro accept a input parameter "barrel" or "endcap" as well, you can do in the prompt line
   for example $ root -l -b -q pngProducer.C("barrel")'. The macro create two png's per chamber the efficiency and the cluster size distributions with the respectively fits. 

4. If you need to get a summary of the fit results for all chambers, you can run the MakeASummary.C macro   just do $ root -l 'MakeASummary.C("barrel", false)' for example,  and you will get a root file in the summary directory with all relevant parameters from the fit. The second bool is the "black-list" parameter, most probably you have to define a black list with problematic rolls in the run and create a txt file "blacklist.txt" with them, if that was the case and you define a black-list with the  help of the experts, and you need a summary only with the filtered rolls(no problematic rolls) so you can put 'true'. 

NOTE:  The blacklist.txt file have to be in the data directory. 

5. The final step of the analysis is a summary with the 7 main parameters. To create the 7 main plots in a special format you need the summary root files ("barrel" and "endcap") and run DrawingOUFlow.C. This macro create 7 png files from 7 TCanvas. Do root -l DrawingOUFlow.C and you get it. 

#Other studies 
 
