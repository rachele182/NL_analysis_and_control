### <font color="green"> <span style="font-size:larger;"> Contents of NL Analysis and Control: </font> </span>

This folder contains all the files used to study, simulate and control the non-linear  system.  
Please refer to the following description:   

- **Simulation Files**:
    - Feedback_dyn.slx simulink file with control of the linearized system (with Input-Output feedback linearization);
    - Feedback_dyn_outLoop.slx simulink file containaing also an outer loop control for the zero dynamics; 
    - Backstepping_control.slx simulink with backstepping control of the rocket + outer loop to stabilize x-axis dynamics;
    - Backstepping_control.slx simulink file with backstepping control complete with estimation of uncertain dynamic parameters i.e mass, inertia and friction coefficients.

- **Matlab Scripts**:
    - init.m: contain the parameters and values to seupt the workspace needed for the simulations; 
    - LIN_analysis.m: script to study the controllability and observability of the system via lienarization method (sufficient but not necessary condition);
    - NL_analysis_contr.m: script to implement filtration method and study local controllability of NL system;
    - NL_analysis_obsv.m: script to implement filtration method and study local observability of NL system;
    - identificability_NL_sys.m: script to study the identificability of friction coefficients;
    - Feedback_LIN.m: script to study if the NL is partially linearizable via dynamic FB;
    - DYN_FL.m: script used to FB-linearize the system: it contains the variables change needed and the control of the linearized system.

- **How to run**:
    1. run matlab script init.m: to setup the work space;
    2. run LIN_analysis.m, NL_analysis_contr.m, NL_analysis_obsv.m: to visualize the results of the analysis with the chosen inputs and outputs of the system;
    3. run identifacibility_NL_sys: to visualize the identificability of the chosen parameters;
    4. run Feddback_LIN.m and DYN_FL.m: to have the result of the complete input-output FB linearization of the system;
    5. open simulink files and run simulation to visualize the results of the control techniques.  

Please note that all the files are written in Matlab and simulink.  
To run the simulations a Simulink version **>=R2020a** is needed.  
