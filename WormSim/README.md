## WormSim c elegans simulation engine from Boyle, Berri and Cohen 2012

See Boyle, Berri and Cohen, [Gait modulation in C. elegans: an integrated neuromechanical model](http://www.frontiersin.org/Computational_Neuroscience/10.3389/fncom.2012.00010/abstract), Front. Comput. Neurosci., 2012.

This file is based on the original readme.txt in the [Supplemental data](http://journal.frontiersin.org/Article/DownloadFile/48074/octet-stream/WormSim1.ZIP/9).

### Instructions for running "WormSim" software on a Linux system.

Installation guide:

1) Clone this project:
    
    git clone https://github.com/OpenSourceBrain/CelegansNeuromechanicalGaitModulation.git

2) Go to "WormSim" directory.

    cd WormSim

3) Set idaInstallDir:

    idaInstallDir=`pwd`/Sundials

4) Go to sundials-2.3.0 directory, configure & make:

    cd sundials-2.3.0
    ./configure CC=g++ --prefix=$idaInstallDir --disable-mpi --disable-fcmix
    make
    make install

5) Look carefully at the resulting output and check for any error messages (there shouldn't be any...). Fix and repeat step 4 if necessary.

6) You should now be able to proceed to running the simulator.

### Usage guide:

1) The whole integrated model is in [worm.cc](https://github.com/OpenSourceBrain/CelegansNeuromechanicalGaitModulation/blob/master/WormSim/Model/worm.cc). Edit this to change any aspects of the model e.g. different environments (water -> intermediate -> agar, with and without objects).

2) To compile the program, make sure you're in the "Model" directory.

3) Enter the command "make".

4) Check the resulting output for errors (some warnings are unfortunately inevitable, but they shouldn't be a problem), fix and re-run "make" if necessary.

5) Enter the command "./program" and wait for it to complete.

6) "program" will generate a file called "simdata.csv", and possibly also "objects.csv" (if objects are being used).

7) Open the appropriate [viewer program](https://github.com/OpenSourceBrain/CelegansNeuromechanicalGaitModulation/blob/master/WormSim/MatlabSupport/WormView.m), found in "WormSim/MatlabSupport/", 
in Matlab. If you don't have access to Matlab, try [Octave](https://www.gnu.org/software/octave/it), or should be possible to reproduce this viewer 
in another language by examining the code and translating it as appropriate.

8) Run the viewer to visualize the model behaviour.

9) Summary:

    cd ../Model
    make
    ./program
    cd ../MatlabSupport
    matlab WormView.m


