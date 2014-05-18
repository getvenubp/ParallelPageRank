Program has been implemented in java language by using MPJ Express
-----------------------------------------
Steps to run from windows(Eclipse)
-----------------------------------------
1. Download MPJ Express(mpj-v0_38) and unpack it.
2. Create a new project and selct mpi.jar from mpj-v0_38\mpi.jar
3. Go to Run -> Run Configurations. Select Environment. 
3. Select new and create new Environment by name MPJ_HOME. Click on Variables. Then Click on Edit Variables. 
	Create new string Substitution, name as MPJ_HOME and values as mpj-v0_38 floder path. Click on ok and select MPJ_HOME from list and click ok.
4. Now Select Argument in Run Configurations.
5. In VM Arguments add, -jar ${MPJ_HOME}\lib\starter.jar
6. In Program Argument add, -np 4 pagerank.input.1000.15 pagerank_output.txt 10 0.001
7. Click on Apply and Run

-----------------------------------------
Steps to run from windows in command line
-----------------------------------------
1. Download MPJ Express and unpack it. 
2. Set MPJ_HOME and PATH environmental variables.
      Right-click My Computer->Properties->Advanced tab->Environment Variables 
      Create new Environment Variable with name as MPJ_HOME and value  as c:\mpj 
      Edit PATH Variable and append the c:\mpj\bin to it. 
3. Compile: javac -cp .;%MPJ_HOME%/lib/mpj.jar PageRankMPI.java
4. Execute: mpjrun.bat -np 4 PageRankMPI pagerank.input.1000.15 pagerank_output.txt 10 0.0001

-----------------------------------------
Steps to run from linux
-----------------------------------------
1. Download MPJ Express and unpack it. 
2. Set MPJ_HOME and PATH environmental variables:
       export MPJ_HOME=/path/to/mpj/
       export PATH=$PATH:$MPJ_HOME/bin 
       (These above two lines can be added to ~/.bashrc)
3. Compile: javac -cp .:$MPJ_HOME/lib/mpj.jar PageRankMPI.java
4. Execute: mpjrun.sh -np 4 PageRankMPI pagerank.input.1000.15 pagerank_output.txt 10 0.0001

	    (mpjrun.sh -np 4 PageRankMPI <Input File Name> <Output File name> <Iteration> <Threshold>)	
