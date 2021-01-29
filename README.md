# Current Propsed Method:
Detailed flow diagram:  
https://docs.google.com/presentation/d/1-b0ZIQYIp3QvBW7MuMyHOe1E7ZB5i-zVTKkM1gTNrMM/edit#slide=id.gb8f5457b7d_0_7

## General Order
1. tree_to_amptools to make post-DSelector trees amptools-ready
2. fullFit_spawn.sh 
3. bootFit_spawn.sh
4. python overlayBins.py 
5. python runPlotEtaPiDeltas.py
6. python diagnoseBS.py
7. python plotIntensities.py

## Example full analysis shell script:
runAnalysis.sh  

Note:  
overlayBins.C uses etapi_plotter program which is not in halld_sim. The code is attached as a folder here that you can insert into halld_sim
