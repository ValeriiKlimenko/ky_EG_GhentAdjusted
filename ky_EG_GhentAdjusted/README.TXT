
The directory "eventGenerator" contains K+Lambda and K+Sigma0
event generators. See eventGenerator/README.TXT for details.

The directory "generalUtils" contains general purpose routines.
It is required to compile the event generator.


This is the event generator for the K^+\Lambda and K^+\Sigma^0
electroproduction in the kinematics region
  0.01 GeV2 < Q2 ~< 12 GeV,
   1.6 GeV  < W  ~< 5 GeV.
The event generator is based on the Jan's Ryckebusch model
(Ghent University)
The output file is in lund format.


To use KLambda of KSigma event generators:
1. install root (https://root.cern.ch/building-root) or type command: use root/6.10.02
2. git clone https://github.com/ValeriiKlimenko/ky_EG_GhentAdjusted
3. cd ky_EG_GhentAdjusted/ky_EG_GhentAdjusted/eventGenerator/
4. chmod +x COMPILE_eg_ky
5. Compile it: 
> ./COMPILE_eg_ky

6. Adjust eg_config_test.txt file. 
 Its format is as follows.
 - first line must contains the string "KLambda" or "KSigma".
 - Next line contains two numbers Q2min and Q2max.
   This is the range in Q2, where the events will be generated.
 - Next line contains two numbers Wmin and Wmax.
   This is the range in W, where the events will be generated.
 - Next line contains the number of events to be generated.
 - Next line must contain the string "dec" or "nodec".
   "dec" means that (e, K+, p, pi-, (gamma)) are the final state
   particles in the output file. Lambda/Sigma0 decay is simulated
   with the phase space approximation.
   "nodec" means that (e, K+, Lambda/Sigma0) are the final state
   particles in the output file.
 - Next line contains the name of the output file.
 - Next line contains the path name to the "data" directory.
   
7. Run it:
>  ./eg_ky eg_config_test.txt
  Parameter eg_config_test.txt is the name of the config. file

Contact person: golovach@jlab.org


