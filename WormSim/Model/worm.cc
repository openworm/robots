/* 
This software is freely available for use in research, teaching and other
non-commercial purposes.  Users have the right to modify, alter, improve,
or enhance the software without limitation, under the condition that you
do not remove or alter the copyright information in this section.

If you publish results which were obtained using the software, or its
source code, please cite the software with a reference to the associated 
publication:

J.H. Boyle, S. Berri and N. Cohen (2012), Gait modulation in C. elegans:
an integrated neuromechanical model, Front. Comput. Neurosci, 6:10 
doi: 10.3389/fncom.2012.00010

Licence information for Sundials IDA can be found in "sundials-2.3.0/LICENCE"
*/

// Required includes
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

using namespace std;

// Simulation paramaters				
#define DURATION  15				//duration of simulation (in seconds)
#define MEDIUM 1				//change in the range 0.0 (water) to 1.0 (agar)
#define OBJECTS 0				//set number of objects (>= 0)
#define LAYOUT 0				//change between 0 (square) 1 (hex) 2 (random)

// Simulator constants
#define NSEG 48
#define NBAR NSEG+1
#define NSEG_MINUS_1 NSEG-1
#define NEQ   3*(NBAR)
#define DELTAT 0.001
#define HALFPI M_PI/2.0

//worm constants
#define NSENSORS 1

// Neural paramaters - Why are these not #defines?
const int N_units = 12;		// Number of neural units
float Hyst = 0.5;		// Neural hysteresis
//const float I_on = 0.675;       // AVB input current (makes the model go)
float I_on = 0;
// General body constants
realtype D = 80e-6;
realtype R[NBAR];
realtype L_seg = 1e-3/NSEG;

// Horizontal element constants
realtype k_PE = (NSEG/24.0)*RCONST(10.0e-3);  	
realtype D_PE = RCONST(0.025)*k_PE;		
realtype AE_PE_ratio = RCONST(20.0);    
realtype k_AE = AE_PE_ratio*k_PE;
realtype D_AE = RCONST(5.0)*AE_PE_ratio*D_PE;	

// Diagonal element constants
realtype k_DE = RCONST(350.0)*k_PE;	
realtype D_DE = RCONST(0.01)*k_DE;	

// Length constant holders
realtype L0_P[NSEG];
realtype L_min[NSEG] ;
realtype L0_P_minus_L_min[NSEG];
realtype L0_D[NSEG];

// Stretch receptor constant holders
realtype SR_shape_compensation[NSEG];

// Muscle time constant
realtype T_muscle = RCONST(0.1);	

// Environment constants
realtype CL_water = RCONST(3.3e-6/(2.0*NBAR));
realtype CN_water = RCONST(5.2e-6/(2.0*NBAR));
realtype CL_agar = RCONST(3.2e-3/(2.0*NBAR));
realtype CN_agar = RCONST(128e-3/(2.0*NBAR));

// Environment variables
realtype K_agar = CN_agar/CL_agar;
realtype CN[NBAR];
realtype CL[NBAR];
realtype ContactForce;
int N_objects = OBJECTS;
int root_N = round(sqrt(N_objects));
realtype Objects[OBJECTS][3];
realtype k_Object = k_PE*5;

// Communication variables  
realtype L_SR[NSEG][2];
realtype I_SR[NSEG][2];

// Neuron and muscle state variables
realtype V_muscle[NSEG][2];
realtype V_neuron[NSEG][2];

// Declare output file
ofstream outfile;

// Prototypes of functions called by IDA (Copied from Sundials examples)
void update_external(realtype timenow);
int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector resval, void *rdata);
void update_neurons(realtype timenow);
void update_muscles(realtype timenow);
void update_SR(realtype timenow);
static int grob(realtype t, N_Vector yy, N_Vector yp, realtype *gout, void *g_data);

// Prototypes of private functions (Copied from Sundials examples)
static void SaveOutput(void *mem, realtype t, N_Vector y);
static void PrintFinalStats(void *mem);
static int check_flag(void *flagvalue, char *funcname, int opt);
double randn(double mu, double sigma);

//prototypes of functions added by Mike for NMI
void excite_neurons(float I_D[],float I_V[], int State[][2],float active); //eventually we should stop using the global neuron variables as well...
void update_SR_Currents(float &active); //stretch receptors

//Sensory's data collection functions
void update_I_Ventral(float I_Array[],int sensorNum);
void update_I_Dorsal(float I_Array[], int sensorNum);
void update_W_Ventral(float W_Array[], int sensorNum);
void update_W_Dorsal(float W_Array[], int sensorNum);
	
//globals for Sensory
//ventral and dorsal sensory currents
float I_ALL_V[NSENSORS][N_units];
float I_ALL_D[NSENSORS][N_units];
//ventral and dorsal sensory weights
float W_ALL_V[NSENSORS][N_units];
float W_ALL_D[NSENSORS][N_units];

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

int main(void)
{
  // IDA variables (Copied from Sundials examples)
  void *mem;
  N_Vector yy, yp, avtol;
  realtype rtol, *yval, *ypval, *atval;
  realtype t0, tout, tret;
  int iout, retval, retvalr;  

  mem = NULL;
  yy = yp = avtol = NULL;
  yval = ypval = atval = NULL;

  // Allocate N-vectors (Copied from Sundials examples)
  yy = N_VNew_Serial(NEQ);
  if(check_flag((void *)yy, "N_VNew_Serial", 0)) return(1);
  yp = N_VNew_Serial(NEQ);
  if(check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);
  avtol = N_VNew_Serial(NEQ);
  if(check_flag((void *)avtol, "N_VNew_Serial", 0)) return(1);

  // Create and initialize  y, y', and absolute tolerance vectors (Copied from Sundials examples)
  yval  = NV_DATA_S(yy);
  ypval = NV_DATA_S(yp);  
  rtol = (MEDIUM < 0.015 ? 0.1 : 1)*RCONST(1.0e-12);
  atval = NV_DATA_S(avtol);

  
  for(int i = 0; i < NBAR; ++i){
	// Initialize body in straight line
	yval[i*3] = i*L_seg;
	yval[i*3+1] = RCONST(0.0);
	yval[i*3+2] = M_PI/RCONST(2.0);

	// Initialize derivative values (Copied from Sundials examples)
 	ypval[i*3] = RCONST(0.0);
  	ypval[i*3+1] = RCONST(0.0);
  	ypval[i*3+2] = RCONST(0.0);

	// Set absolute tolerances for solver (Copied from Sundials examples)
	// Tolerance must be set lower when simulating in water, due to lower drag coefficients
	atval[i*3] = (MEDIUM < 0.015 ? 0.1 : 1)*RCONST(1.0e-9);
	atval[i*3+1] = (MEDIUM < 0.015 ? 0.1 : 1)*RCONST(1.0e-9);
	atval[i*3+2] = (MEDIUM < 0.015 ? 0.1 : 1)*RCONST(1.0e-5);

  }
  
  // Initialize model variables
  for(int i = 0; i < NSEG; ++i){
	V_muscle[i][0] = 0.0;
	V_muscle[i][1] = 0.0;
	V_neuron[i][0] = 0.0;
	V_neuron[i][1] = 0.0;
	I_SR[i][0] = 0.0;
	I_SR[i][1] = 0.0;
  }  
  
  // Set local body radius values based on elliptical approximation
  for(int i = 0; i < NBAR; ++i){	
	R[i] = D/2.0*fabs(sin(acos((i-NSEG/2.0)/(NSEG/2.0 + 0.2))));		
  } 

  // Set stretch receptor weightings that compensate for the elliptical shape,
  // giving approximately the same SR response to segment mending angle
  for(int i = 0; i < NSEG; ++i){
	SR_shape_compensation[i] = D/(R[i] + R[i+1]);  
  } 
 
  // Set muscle constants (rest length, minimum length etc) accounting 
  // for length differences due to the elliptical body shape
  for(int i = 0; i < NSEG; ++i){
	float scale = 0.65*((R[i] + R[i+1])/D);
	L0_P[i] = sqrt(pow(L_seg,2) + pow((R[i] - R[i+1]),2));
	L_min[i] = (1.0-scale)*L0_P[i];
	L0_P_minus_L_min[i] = L0_P[i] - L_min[i];
	L0_D[i] = sqrt(pow(L_seg,2) + pow((R[i] + R[i+1]),2));
  }

  // Set drag constants according to medium
  for(int i = 0; i < NBAR; ++i){
  	CL[i] = (CL_agar - CL_water)*MEDIUM + CL_water;
  	CN[i] = (CN_agar - CN_water)*MEDIUM + CN_water;
  }
		

  // Place objects in environment (if using them)  
  if(N_objects > 0){
	if(LAYOUT == 0){
  		//Square post array ala Park et al.
		//Adjust these parameters to modify configuration
  		float post_radius = 0.1e-3;
  		float post_spacing = 0.4e-3;
  		ofstream object_outfile;
  		object_outfile.open("objects.csv");
  		for(int i = 0; i < root_N; ++i){
			for(int j = 0; j < root_N; ++j){
				Objects[root_N*i + j][0] = post_spacing*j - 0.75*(root_N-1)*post_spacing;
				Objects[root_N*i + j][1] = post_spacing*i - 0.25*(root_N-1)*post_spacing;
				Objects[root_N*i + j][2] = post_radius;
				object_outfile << Objects[root_N*i + j][0] << ", " << Objects[root_N*i + j][1] << ", " << Objects[root_N*i + j][2] << "\n";
			}
  		}
  		object_outfile.close();	
	}
	else if(LAYOUT == 1){
		//Hexagonal post array ala Lockery et al.
		//Adjust these parameters to modify configuration
  		float post_radius = 0.1e-3;
  		float min_gap = 0.2e-3;
  		float horizontal_post_spacing = min_gap + 2*post_radius;
  		float vertical_post_spacing = sqrt(pow(horizontal_post_spacing,2) - pow(0.5*horizontal_post_spacing,2));
  		ofstream object_outfile;
  		object_outfile.open("objects.csv");
  		for(int i = 0; i < root_N; ++i){
			for(int j = 0; j < root_N; ++j){
				Objects[root_N*i + j][0] = horizontal_post_spacing*j - 0.85*(root_N-1)*horizontal_post_spacing + 0.5*horizontal_post_spacing*(i%2 == 0);
				Objects[root_N*i + j][1] = vertical_post_spacing*i - 0.5*(root_N-1)*vertical_post_spacing;
				Objects[root_N*i + j][2] = post_radius;
				object_outfile << Objects[root_N*i + j][0] << ", " << Objects[root_N*i + j][1] << ", " << Objects[root_N*i + j][2] << "\n";
			}
  		}
  		object_outfile.close();
	}
	else if(LAYOUT == 2){
  		//Create random layout of objects
		//Adjust these parameters to modify configuration
  		float X_lim = 1.5e-3;
  		float Y_lim = 1.5e-3;
  		float min_gap = 0.08e-3;
  		float R_min = 0.04e-3;
  		float R_max = 0.25e-3;
  		bool fits;
  		ofstream object_outfile;
  		object_outfile.open("objects.csv");
  		for(int i = 0; i < N_objects; ++i){
			do{
				fits = true;
				Objects[i][0] = -1.5*X_lim + (2*X_lim)*(rand()/(RAND_MAX*1.0));
				Objects[i][1] = -Y_lim + (2*Y_lim)*(rand()/(RAND_MAX*1.0));
				Objects[i][2] = R_min + (R_max - R_min)*(rand()/(RAND_MAX*1.0));
	
				for (int j = 0; j < i; ++j){
					float dist = sqrt(pow((Objects[i][0] - Objects[j][0]),2) + pow((Objects[i][1] - Objects[j][1]),2));
					if(dist < (Objects[i][2] + Objects[j][2] + min_gap)){
						fits = false;
					}
				}
				for(int j = 0; j < NBAR; ++j){
					float dist = sqrt(pow((Objects[i][0] - yval[j*3]),2) + pow((Objects[i][1] - yval[j*3+1]),2));
					if(dist < (Objects[i][2] + 0.05e-3)){
						fits = false;
					}
				}
			}
			while(fits == false);
			object_outfile << Objects[i][0] << ", " << Objects[i][1] << ", " << Objects[i][2] << "\n";
  		}
  		object_outfile.close();
	} 
  	cerr << "\nObjects set up";
  }

  // Integration start time 
  t0 = RCONST(0.0);
  
  // Call IDACreate and IDAMalloc to initialize IDA memory (Copied from Sundials examples)
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);

  retval = IDAMalloc(mem, resrob, t0, yy, yp, IDA_SV, rtol, avtol);	
  if(check_flag(&retval, "IDAMalloc", 1)) return(1);
  
  // Free avtol (Copied from Sundials examples) 
  N_VDestroy_Serial(avtol);

  // Call IDADense and set up the linear solver (Copied from Sundials examples) 
  retval = IDADense(mem, NEQ);
  if(check_flag(&retval, "IDADense", 1)) return(1);

  // Open output file  
  outfile.open("simdata.csv");  
  
  // Integrator inputs
  iout = 0; 
  tout = DELTAT;

  // Save current state
  SaveOutput(mem,0,yy);

  //these simulations take a WHILE, so let's add a timing output
  double percDone = tout/DURATION*100;
  int nPers = 1;

  // Loop through model for specified amount of time
  while(tout <= DURATION) {
   
    	// Call residual function (Copied from Sundials examples)
	// to update physical model (Sundials takes multiple steps)
    	retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);    	
    
    	// Call stretch receptor update function
    	update_SR(tout);

    	// Call neural model update function
    	update_neurons(tout);

    	//Call muscle model update function
    	update_muscles(tout);

    	// Save current state
    	SaveOutput(mem,tret,yy);

	// Check integration went ok (Copied from Sundials examples) 
    	if(check_flag(&retval, "IDASolve", 1)) return(1);
    
	// Prepare to go to next step
    	if (retval == IDA_SUCCESS) {
      		iout++;
      		tout += DELTAT;
    	}
	if(percDone > nPers*10)
	{
		nPers++;
		cout << percDone << " Percent done." << endl;
	}
	percDone = tout/DURATION*100;
  }

  // (Copied from Sundials examples) 
  PrintFinalStats(mem);

  // Free memory (Copied from Sundials examples) 
  IDAFree(&mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);

  // Close output file
  outfile.close();

  return(0);
  
}

/*
 *--------------------------------------------------------------------
 * Model Functions
 *--------------------------------------------------------------------
 */



//Sensory's update helper function
void update_Sensory_Data(float inAr[], float destAr[][N_units], int sensorNum)
{
	for(int i = 0; i < N_units; i++)
	{
		destAr[sensorNum][i] = inAr[i];
	}
}

//Sensory's data collection functions
void update_I_Ventral(float I_Array[],int sensorNum)
{
	update_Sensory_Data(I_Array,I_ALL_V,sensorNum);
}
void update_I_Dorsal(float I_Array[], int sensorNum)
{
	update_Sensory_Data(I_Array, I_ALL_D,sensorNum);
}
void update_W_Ventral(float W_Array[], int sensorNum)
{
	update_Sensory_Data(W_Array, W_ALL_V, sensorNum);
}
void update_W_Dorsal(float W_Array[], int sensorNum)
{
	update_Sensory_Data(W_Array, W_ALL_D,sensorNum);
}

  // Neural circuit function.  This function takes in currents from all Sensors and computes the weighted sum I_D and I_V neuron input currents
  void update_neurons(realtype timenow){
	  
	// Neural state variables
   	static int State[N_units][2];
	
	// If this is the first time update_neurons is called, initialize with all neurons on one side ON - we won't actually be using this once we incorporate Sensory
   	static bool initialized = false;     	
	if(!initialized){
		for(int i = 0; i < N_units; ++i){
			State[i][0] = 1;
			State[i][1] = 0;
		}
		initialized = true;
	}
   
	 

   	// Variables for total input current to each B-class motorneuron
   	float I_D[N_units] = {0};
   	float I_V[N_units] = {0}; //reset to zero, or else the summation will cause serious issues below.    

   	// Current bias to compensate for the fact that neural inhibition only goes one way
   	float I_bias = 0.5;

	// Combine AVB current, stretch receptor current, neural inhibition and bias
	  /*
	  	Mike's note: This "NSENSORS" stuff came about due to a misunderstanding I had.  I still need to go through and fix this part,
		but it will require significant modification of the code.
	  */
	  for(int j = 0; j < NSENSORS; j++)
	  {
			for(int i = 0; i < N_units; ++i){
				I_D[i] = I_D[i] + I_on + W_ALL_D[j][i]*I_ALL_D[j][i];
				I_V[i] = I_V[i] + (I_bias - State[i][j]) + I_on + W_ALL_V[j][i]*I_ALL_V[j][i];

			}
	  }
	  
	  
	  /*
	  	This code below is the *very* beginning of the command neuron circuit.  It enables us to move forward and backward.
		
		The current represents either AVA or AVB, and the variable act is a poweful variable - it represents a muscle voltage scaling factor
		that may allow us to modulate response to touch magnitude, and its sign governs which way we move.  
		
		TODO: This section requires a lot of cleanup
	  */
	  float act = -1;
	  I_on = 0.675;
	  if(timenow < 4)
	  {
		  act = 1;
	  }
	  else if(timenow >= 4 && timenow < 6)
		  I_on = 0;
	  else if(timenow >=6)
	  {
		  act = -1;
	  }
	  
	//now excite
	excite_neurons(I_D,I_V, State,act);
	 
   
  }

/*
	This function takes the summed I_D and I_V and activates neurons.
	
	Active serves two roles.  If it is negative, it reverses direction.  Also, it scales V_neuron to modulate response strength
*/

void excite_neurons(float I_D[],float I_V[], int State[][2], float active)
{
	// Set up neuromuscular junctions   
   	static float NMJ_weight[NSEG];
 
	//initialization depends on direction
	if(active >0)
	{
		for(int i = 0; i < NSEG; ++i){
			NMJ_weight[i] =  0.7*(1.0 - i * 0.6/NSEG);	// Decreasing gradient in NMJ strength / muscle efficacy
		}
		NMJ_weight[0] /= 1.5;				// Helps to prevent excessive bending of head
	}
	else{
		for(int i = 0; i <NSEG; ++i){
			NMJ_weight[i] =  0.7*(1.0 - (NSEG-i-1) * 0.6/NSEG);	// Decreasing gradient in NMJ strength / muscle efficacy
		}
		NMJ_weight[NSEG-1] /= 1.5;				// Helps to prevent excessive bending of head
	}
	
	//active gets modified in this function - if it is negative, it becomes positive.
	update_SR_Currents(active);

	   	// Update state for each bistable B-class neuron
   	for(int i = 0; i < N_units; ++i){
   		if(I_D[i] > (0.5 + Hyst/2.0 - Hyst*State[i][0])){
			State[i][0] = 1;}
  		else{
			State[i][0] = 0;}

   		if(I_V[i] > (0.5 + Hyst/2.0 - Hyst*State[i][1])){
			State[i][1] = 1;}
   		else{
			State[i][1] = 0;}
   	}
	

   	// Compute effective input to each muscle including B-class excitation and contralateral D-class inhibition, scaled by active to modulate response magnitude 	
   	for(int i = 0; i < NSEG; ++i){		
		V_neuron[i][0] = active*(NMJ_weight[i]*State[(int)(i*N_units/NSEG)][0] - NMJ_weight[i]*State[(int)(i*N_units/NSEG)][1]);
		V_neuron[i][1] = active*(NMJ_weight[i]*State[(int)(i*N_units/NSEG)][1] - NMJ_weight[i]*State[(int)(i*N_units/NSEG)][0]);	
   	}
	
}

//Stretch Receptor Functions

// Update the stretch receptors (for each segment). These
// are weighted and combined as input to the neural units in function "update_neurons"
void update_SR(realtype timenow){
for(int i = 0; i < NSEG; ++i){	
	// Bilinear SR function on one side to compensate for asymmetry and help worm go straight
	I_SR[i][0] = SR_shape_compensation[i]*((L_SR[i][0] - L0_P[i])/L0_P[i]*((L_SR[i][0] > L_seg) ? 0.8:1.2));
	I_SR[i][1] = SR_shape_compensation[i]*((L_SR[i][1] - L0_P[i])/L0_P[i]);
}

}

/*
	Updates the currents from the stretch receptors, then sends them off to Sensory
*/
void update_SR_Currents(float &active)
{
		  
   	// Stretch receptor variables   	
   	static float I_SR_D[N_units];
   	static float I_SR_V[N_units];
   	static float SR_weight[N_units];

   	int N_SR = 6;	//This refers to the number of UNITS (not segments) that each unit receives feedback from (thus 1 means just local feedback)   
   	int N_seg_per_unit = (int)(NSEG/N_units);   
	
	 
   	// SR_weight is a global weighting for each unit, used to get the compensate for curvature gradient induced by the NMJ gradient above
   
	//what you do depends on which direction, forward or backward, you are going.
	if(active > 0)
	{
		for(int i = 0; i < N_units; ++i){	
			SR_weight[i] = 0.65*(0.4 + 0.08*i)*(N_units/12.0)*(2.0/N_seg_per_unit);	
		}
	}
	else
	{
		for(int i = 0; i < N_units; ++i){	
			SR_weight[i] = 0.65*(0.4 + 0.08*(N_units-i-1))*(N_units/12.0)*(2.0/N_seg_per_unit);	
		}
	}
  
	//once again, there is a subtle difference between forward and backward movement
	if(active < 0)
	{
		//now reset active for use in functions above
		active = -active;
		for(int i = N_SR; i <= N_units; ++i){
			I_SR_D[i] = 0.0;
			I_SR_V[i] = 0.0;
			for(int j = 0; j < N_SR; ++j){
				I_SR_D[i] += I_SR[(i-j)*N_seg_per_unit][0] + (N_seg_per_unit >= 2)*I_SR[(i-j)*N_seg_per_unit + 1][0] + (N_seg_per_unit >= 3)*I_SR[(i-j)*N_seg_per_unit + 2][0]  + (N_seg_per_unit >= 4)*I_SR[(i-j)*N_seg_per_unit + 3][0];	
				I_SR_V[i] += I_SR[(i-j)*N_seg_per_unit][1] + (N_seg_per_unit >= 2)*I_SR[(i-j)*N_seg_per_unit + 1][1] + (N_seg_per_unit >= 3)*I_SR[(i-j)*N_seg_per_unit + 2][1]  + (N_seg_per_unit >= 4)*I_SR[(i-j)*N_seg_per_unit + 3][1];				
			}		
		}

		// For units near the tail, fewer segments contribute (because the body ends)
		int tmp_N_SR = N_SR;
		for(int i = (N_SR - 1); i >= 0; --i){			
			tmp_N_SR --;
			I_SR_D[i] = 0.0;
			I_SR_V[i] = 0.0;
			for(int j = 0; j < tmp_N_SR; ++j){
				I_SR_D[i] += I_SR[(i-j)*N_seg_per_unit][0] + (N_seg_per_unit >= 2)*I_SR[(i-j)*N_seg_per_unit + 1][0] + (N_seg_per_unit >= 3)*I_SR[(i-j)*N_seg_per_unit + 2][0]  + (N_seg_per_unit >= 4)*I_SR[(i-j)*N_seg_per_unit + 3][0];	
				I_SR_V[i] += I_SR[(i-j)*N_seg_per_unit][1] + (N_seg_per_unit >= 2)*I_SR[(i-j)*N_seg_per_unit + 1][1] + (N_seg_per_unit >= 3)*I_SR[(i-j)*N_seg_per_unit + 2][1]  + (N_seg_per_unit >= 4)*I_SR[(i-j)*N_seg_per_unit + 3][1];
			}			
		}	

		// Compensate for the anterior segments with shorter processes   	
		for(int i = (N_SR - 1); i >= 0; --i){ 
			I_SR_D[i] *= sqrt((N_SR/(i+1)));
			I_SR_V[i] *= sqrt((N_SR/(i+1)));	
		}
	}
	else
	{
		// Add up stretch receptor contributions from all body segments in receptive field for each neural unit 
		for(int i = 0; i <= N_units - N_SR; ++i){
			I_SR_D[i] = 0.0;
			I_SR_V[i] = 0.0;
			for(int j = 0; j < N_SR; ++j){
				I_SR_D[i] += I_SR[(i+j)*N_seg_per_unit][0] + (N_seg_per_unit >= 2)*I_SR[(i+j)*N_seg_per_unit + 1][0] + (N_seg_per_unit >= 3)*I_SR[(i+j)*N_seg_per_unit + 2][0]  + (N_seg_per_unit >= 4)*I_SR[(i+j)*N_seg_per_unit + 3][0];	
				I_SR_V[i] += I_SR[(i+j)*N_seg_per_unit][1] + (N_seg_per_unit >= 2)*I_SR[(i+j)*N_seg_per_unit + 1][1] + (N_seg_per_unit >= 3)*I_SR[(i+j)*N_seg_per_unit + 2][1]  + (N_seg_per_unit >= 4)*I_SR[(i+j)*N_seg_per_unit + 3][1];				
			}		
		}

		// For units near the tail, fewer segments contribute (because the body ends)
		int tmp_N_SR = N_SR;
		for(int i = (N_units - N_SR + 1); i < N_units; ++i){			
			tmp_N_SR --;
			I_SR_D[i] = 0.0;
			I_SR_V[i] = 0.0;
			for(int j = 0; j < tmp_N_SR; ++j){
				I_SR_D[i] += I_SR[(i+j)*N_seg_per_unit][0] + (N_seg_per_unit >= 2)*I_SR[(i+j)*N_seg_per_unit + 1][0] + (N_seg_per_unit >= 3)*I_SR[(i+j)*N_seg_per_unit + 2][0]  + (N_seg_per_unit >= 4)*I_SR[(i+j)*N_seg_per_unit + 3][0];	
				I_SR_V[i] += I_SR[(i+j)*N_seg_per_unit][1] + (N_seg_per_unit >= 2)*I_SR[(i+j)*N_seg_per_unit + 1][1] + (N_seg_per_unit >= 3)*I_SR[(i+j)*N_seg_per_unit + 2][1]  + (N_seg_per_unit >= 4)*I_SR[(i+j)*N_seg_per_unit + 3][1];
			}			
		}	

		// Compensate for the posterior segments with shorter processes   	
		for(int i = (N_units - N_SR + 1); i < N_units; ++i){ 
			I_SR_D[i] *= sqrt(-(N_SR/(i-N_units)));
			I_SR_V[i] *= sqrt(-(N_SR/(i-N_units)));	
		}
	}
	//update Sensory
	update_I_Ventral(I_SR_V,0);
	update_I_Dorsal(I_SR_D,0);
	update_W_Ventral(SR_weight,0);
	update_W_Dorsal(SR_weight,0);
}

  // Update the simple muscle "model" (electronic)
  void update_muscles(realtype timenow){
  	//Muscle transfer function is just a simple LPF
  	for(int i = 0; i < NSEG; ++i){
		for(int j = 0; j < 2; ++j){
			realtype dV = (V_neuron[i][j] - V_muscle[i][j])/T_muscle;
			V_muscle[i][j] += dV*DELTAT;		
		}
  	}
  }

  // System residual function which implements physical model (Based on Sundials examples) 
  int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *rdata)
  {
  	// Import data from vectors
  	realtype *yval, *ypval, *rval;
  	yval = NV_DATA_S(yy); 
  	ypval = NV_DATA_S(yp); 
  	rval = NV_DATA_S(rr);  

  	//Declare variables
  	realtype CoM[NBAR][3];
  	realtype V_CoM[NBAR][3];
  	realtype term[NBAR][2][2];  		// Nseg, d/v, x/y
  	realtype V_term[NBAR][2][2];
  	realtype dy,dx,dVy,dVx,F_even,F_odd;
  	realtype F_term[NBAR][2][2];
  	realtype F_term_rotated[NBAR][2][2];
  	realtype V_CoM_rotated[NBAR][3];

  	realtype L[NSEG][2];				// Nseg, d/v
  	realtype Dir[NSEG][2][2];			// Nseg, d/v, x/y
  	realtype S[NSEG][2];
  	realtype L_D[NSEG][2];			// Nseg, \,/   <- these are the angles of the diagonals
  	realtype Dir_D[NSEG][2][2];  			// Nseg, \,/ , x/y
  	realtype S_D[NSEG][2];

  	realtype L0_AE, T, F_AE, F_PE, F_PD;
  	realtype F_H[NSEG][2];  
  	realtype F_D[NSEG][2];

  	realtype L_EXT;
  	realtype F_EXT[2];

  	realtype F_object[NBAR][2][2];
	// Initialize all object forces to zero incase objects are not being used
	for(int i = 0; i < NBAR; ++i){
		for(int j = 0; j < 2; ++j){
			for(int k = 0; k < 2; ++k){
				F_object[i][j][k] = 0;
			}
		}
	}
  	
  	for(int i = 0; i < NBAR; ++i){
		// Extract CoM of each solid rod from vectors
		int three_i = i*3;	
		CoM[i][0] = yval[three_i];
		CoM[i][1] = yval[three_i + 1];
		CoM[i][2] = yval[three_i + 2];

		// Calculate positions of D/V points based on CoM, angle and radius
		dx = R[i]*cos(CoM[i][2]);
		dy = R[i]*sin(CoM[i][2]);

		term[i][0][0] = CoM[i][0] + dx;
		term[i][0][1] = CoM[i][1] + dy;
		term[i][1][0] = CoM[i][0] - dx;
		term[i][1][1] = CoM[i][1] - dy;

		// Extract CoM velocities of each solid rod from vectors
		V_CoM[i][0] = ypval[three_i];
		V_CoM[i][1] = ypval[three_i + 1];
		V_CoM[i][2] = ypval[three_i + 2];

		// Calculate velocity of D/V points based on CoM velocity, rate of rotation and radius
		realtype V_arm = R[i]*V_CoM[i][2];		
		dVx = V_arm*cos(CoM[i][2] + HALFPI);  		
		dVy = V_arm*sin(CoM[i][2] + HALFPI);  		

		V_term[i][0][0] = V_CoM[i][0] + dVx;
		V_term[i][0][1] = V_CoM[i][1] + dVy;
		V_term[i][1][0] = V_CoM[i][0] - dVx;
		V_term[i][1][1] = V_CoM[i][1] - dVy;		
  	}	

  
  	// Get Horizontal/Diagonal element lengths and lengthening/shortening velocities
  	for(int i = 0; i < NSEG; ++i){

		// Strange format for efficiency
		int iplus1 = i+1;
	
		Dir[i][0][0] = (term[iplus1][0][0] - term[i][0][0]);
		Dir[i][0][1] = (term[iplus1][0][1] - term[i][0][1]);
		L[i][0] = sqrt(pow(Dir[i][0][0],2.0) + pow(Dir[i][0][1],2.0));
		Dir[i][0][0] /= L[i][0];
		Dir[i][0][1] /= L[i][0];
		S[i][0] = (V_term[iplus1][0][0] - V_term[i][0][0])*Dir[i][0][0] + (V_term[iplus1][0][1] - V_term[i][0][1])*Dir[i][0][1];

		Dir[i][1][0] =  (term[iplus1][1][0] - term[i][1][0]);
		Dir[i][1][1] =  (term[iplus1][1][1] - term[i][1][1]);
		L[i][1] = sqrt(pow(Dir[i][1][0],2.0) + pow(Dir[i][1][1],2.0));
		Dir[i][1][0] /= L[i][1];
		Dir[i][1][1] /= L[i][1];
		S[i][1] = (V_term[iplus1][1][0] - V_term[i][1][0])*Dir[i][1][0] + (V_term[iplus1][1][1] - V_term[i][1][1])*Dir[i][1][1];

		Dir_D[i][0][0] =  (term[iplus1][1][0] - term[i][0][0]);
		Dir_D[i][0][1] =  (term[iplus1][1][1] - term[i][0][1]);
		L_D[i][0] = sqrt(pow(Dir_D[i][0][0],2.0) + pow(Dir_D[i][0][1],2.0));
		Dir_D[i][0][0] /= L_D[i][0];
		Dir_D[i][0][1] /= L_D[i][0];
		S_D[i][0] = (V_term[iplus1][1][0] - V_term[i][0][0])*Dir_D[i][0][0] + (V_term[iplus1][1][1] - V_term[i][0][1])*Dir_D[i][0][1];

		Dir_D[i][1][0] =  (term[iplus1][0][0] - term[i][1][0]);
		Dir_D[i][1][1] =  (term[iplus1][0][1] - term[i][1][1]);
		L_D[i][1] = sqrt(pow(Dir_D[i][1][0],2.0) + pow(Dir_D[i][1][1],2.0));
		Dir_D[i][1][0] /= L_D[i][1];
		Dir_D[i][1][1] /= L_D[i][1];
		S_D[i][1] = (V_term[iplus1][0][0] - V_term[i][1][0])*Dir_D[i][1][0] + (V_term[iplus1][0][1] - V_term[i][1][1])*Dir_D[i][1][1];  

  		// Calculate force contributions on each D/V point

  		//Dorsal forces due to horizontal elements
  		L0_AE = L0_P[i] - fmax(V_muscle[i][0],0)*(L0_P_minus_L_min[i]);  	
  		
  		F_AE = k_AE*fmax(V_muscle[i][0],0)*(L0_AE - L[i][0]);
  		F_PE = k_PE*((L0_P[i] - L[i][0]) + ((L[i][0]-L0_P[i]) > RCONST(0.0))*pow(RCONST(2.0)*(L[i][0]-L0_P[i]),4));
  		F_PD = (D_PE + fmax(V_muscle[i][0],0)*D_AE)*S[i][0];

  		F_H[i][0] = F_PE + F_AE - F_PD;
	  
  		//Ventral forces due to horizontal elements
  		L0_AE = L0_P[i] - fmax(V_muscle[i][1],0)*(L0_P_minus_L_min[i]);  	
  		
  		F_AE = k_AE*fmax(V_muscle[i][1],0)*(L0_AE - L[i][1]);
  		F_PE = k_PE*((L0_P[i] - L[i][1]) + ((L[i][1]-L0_P[i]) > RCONST(0.0))*pow(RCONST(2.0)*(L[i][1]-L0_P[i]),4));
  		F_PD = (D_PE + fmax(V_muscle[i][1],0)*D_AE)*S[i][1];

  		F_H[i][1] = F_PE + F_AE - F_PD;
	
  		//Diagonal forces due to diagonal elements
  		F_D[i][0] = (L0_D[i] - L_D[i][0])*k_DE - D_DE*S_D[i][0];
  		F_D[i][1] = (L0_D[i] - L_D[i][1])*k_DE - D_DE*S_D[i][1];
  	}

  	// If using objects, check for object collisions and calculate associated forces
	if(N_objects > 0){
  		realtype P_x,P_y,Distance,magF,D_scale,magF_P1,magF_P2;
  		ContactForce = 0;
  		for(int i = 0; i < NBAR; ++i){
			for(int j = 0; j < 2; ++j){
				// First ensure they contain zeros
				F_object[i][j][0] = 0;
				F_object[i][j][1] = 0;
				P_x = term[i][j][0];
				P_y = term[i][j][1];

				// Now check proximity to each object
				for(int k = 0; k < N_objects; ++k){
					if((P_x<(Objects[k][0]+Objects[k][2]))&&(P_x>(Objects[k][0]-Objects[k][2]))&&(P_y<(Objects[k][1]+Objects[k][2]))&&(P_y>(Objects[k][1]-Objects[k][2]))){

						//This means the point is within the bounding box of the object, so now we must compute the force (if any)
						dx = P_x - Objects[k][0];
						dy = P_y - Objects[k][1];
						Distance = sqrt(pow(dx,2) + pow(dy,2));
						D_scale = 0.01*Objects[k][2];

						if(Distance < Objects[k][2]){
							magF = k_Object*(Objects[k][2] - Distance) + D_scale*k_Object*(pow((Objects[k][2] - Distance)/D_scale,2));
							F_object[i][j][0] += (dx/Distance)*magF;
							F_object[i][j][1] += (dy/Distance)*magF;
							ContactForce += magF;
						}
					}			
				}
			}
  		}
	}

  	// Add up force contributions for each D/V point
  	F_term[0][0][0] = -F_H[0][0]*Dir[0][0][0] - F_D[0][0]*Dir_D[0][0][0] + F_object[0][0][0];
  	F_term[0][0][1] = -F_H[0][0]*Dir[0][0][1] - F_D[0][0]*Dir_D[0][0][1] + F_object[0][0][1];

  	F_term[0][1][0] = -F_H[0][1]*Dir[0][1][0] - F_D[0][1]*Dir_D[0][1][0] + F_object[0][1][0];
  	F_term[0][1][1] = -F_H[0][1]*Dir[0][1][1] - F_D[0][1]*Dir_D[0][1][1] + F_object[0][1][1];

  	for(int i = 1; i < NSEG; ++i){
		int i_minus_1 = i-1;

		F_term[i][0][0] = F_H[i_minus_1][0]*Dir[i_minus_1][0][0] - F_H[i][0]*Dir[i][0][0] + F_D[i_minus_1][1]*Dir_D[i_minus_1][1][0] - F_D[i][0]*Dir_D[i][0][0] + F_object[i][0][0];
		F_term[i][0][1] = F_H[i_minus_1][0]*Dir[i_minus_1][0][1] - F_H[i][0]*Dir[i][0][1] + F_D[i_minus_1][1]*Dir_D[i_minus_1][1][1] - F_D[i][0]*Dir_D[i][0][1] + F_object[i][0][1];

		F_term[i][1][0] = F_H[i_minus_1][1]*Dir[i_minus_1][1][0] - F_H[i][1]*Dir[i][1][0] + F_D[i_minus_1][0]*Dir_D[i_minus_1][0][0] - F_D[i][1]*Dir_D[i][1][0] + F_object[i][1][0];
		F_term[i][1][1] = F_H[i_minus_1][1]*Dir[i_minus_1][1][1] - F_H[i][1]*Dir[i][1][1] + F_D[i_minus_1][0]*Dir_D[i_minus_1][0][1] - F_D[i][1]*Dir_D[i][1][1] + F_object[i][1][1];
  	}

  	F_term[NSEG][0][0] = F_H[NSEG_MINUS_1][0]*Dir[NSEG_MINUS_1][0][0] + F_D[NSEG_MINUS_1][1]*Dir_D[NSEG_MINUS_1][1][0] + F_object[NSEG][0][0];
  	F_term[NSEG][0][1] = F_H[NSEG_MINUS_1][0]*Dir[NSEG_MINUS_1][0][1] + F_D[NSEG_MINUS_1][1]*Dir_D[NSEG_MINUS_1][1][1] + F_object[NSEG][0][1];

  	F_term[NSEG][1][0] = F_H[NSEG_MINUS_1][1]*Dir[NSEG_MINUS_1][1][0] + F_D[NSEG_MINUS_1][0]*Dir_D[NSEG_MINUS_1][0][0] + F_object[NSEG][1][0];
  	F_term[NSEG][1][1] = F_H[NSEG_MINUS_1][1]*Dir[NSEG_MINUS_1][1][1] + F_D[NSEG_MINUS_1][0]*Dir_D[NSEG_MINUS_1][0][1] + F_object[NSEG][1][1];
  
  	// Convert net forces on D/V points to force and torque	acting on rod CoM
  	for(int i = 0; i < NBAR; ++i){
		realtype cos_thi = cos(CoM[i][2]);
		realtype sin_thi = sin(CoM[i][2]);
		for(int j = 0; j < 2; ++j){			
			F_term_rotated[i][j][0] = F_term[i][j][0]*cos_thi + F_term[i][j][1]*sin_thi;	// This is Fperp
			F_term_rotated[i][j][1] = F_term[i][j][0]*sin_thi - F_term[i][j][1]*cos_thi;    // THis is Fparallel
		}

		V_CoM_rotated[i][0] = (F_term_rotated[i][0][0] + F_term_rotated[i][1][0])/CN[i];

		F_even = (F_term_rotated[i][0][1] + F_term_rotated[i][1][1]);	//Took out the /2
		F_odd = (F_term_rotated[i][1][1] - F_term_rotated[i][0][1])/RCONST(2.0);	

		V_CoM_rotated[i][1] = (F_even)/CL[i];				//Allowing me to take out *2
		V_CoM[i][2] = (F_odd/CL[i])/(M_PI*2.0*R[i]);

		V_CoM[i][0] = V_CoM_rotated[i][0]*cos_thi + V_CoM_rotated[i][1]*sin_thi;
		V_CoM[i][1] = V_CoM_rotated[i][0]*sin_thi - V_CoM_rotated[i][1]*cos_thi;
	
		int three_i = i*3;

		rval[three_i] = V_CoM[i][0] - ypval[three_i];
		rval[three_i+1] = V_CoM[i][1] - ypval[three_i+1];
		rval[three_i+2] = V_CoM[i][2] - ypval[three_i+2];
  	}

  	// Store old lengths for Stretch Receptors 
  	for(int i = 0; i < NSEG; ++i){
		L_SR[i][0] = L[i][0];
		L_SR[i][1] = L[i][1];
  	}

  	return(0);
  }


/*
 *--------------------------------------------------------------------
 * Private functions
 *--------------------------------------------------------------------
 */


  // Save output to csv file for visualization
  static void SaveOutput(void *mem, realtype t, N_Vector y){
  
  	float FPS = 25;	// Frame rate at which output should be logged
  	static int tCount = 1.0/(DELTAT*FPS);
  	realtype *yval;
  	int retval, kused;
  	long int nst;
  	realtype hused;

  	yval  = NV_DATA_S(y);
  
  	// Save data to file
  	if(tCount == 1.0/(DELTAT*FPS)){
  
  		outfile << t;	
  		for(int i = 0; i < NBAR; ++i){
			outfile << ", " << yval[i*3] << ", " << yval[i*3+1] << ", " << yval[i*3+2];
  		}
  		outfile << "\n";
		tCount = 0;
  	}

  	++tCount;
  	return;  
  }

  double randn(double mu, double sigma) {
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double polar, rsquared, var1, var2;
	
	//	If no deviate has been stored, the polar Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {
		
		//	choose pairs of uniformly distributed deviates, discarding those 
		//	that don't fall within the unit circle
		do {
			var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0 || rsquared == 0.0);
		
		//	calculate polar tranformation for each deviate
		polar=sqrt(-2.0*log(rsquared)/rsquared);
		
		//	store first deviate and set flag
		storedDeviate=var1*polar;
		deviateAvailable=true;
		
		//	return second deviate
		return var2*polar*sigma + mu;
	}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
  }



  // Print final integrator statistics (Copied from Sundials examples)
  static void PrintFinalStats(void *mem){
  	int retval;
  	long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;

  	retval = IDAGetNumSteps(mem, &nst);
  	check_flag(&retval, "IDAGetNumSteps", 1);
  	retval = IDAGetNumResEvals(mem, &nre);
  	check_flag(&retval, "IDAGetNumResEvals", 1);
  	retval = IDADenseGetNumJacEvals(mem, &nje);
  	check_flag(&retval, "IDADenseGetNumJacEvals", 1);
  	retval = IDAGetNumNonlinSolvIters(mem, &nni);
  	check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  	retval = IDAGetNumErrTestFails(mem, &netf);
  	check_flag(&retval, "IDAGetNumErrTestFails", 1);
  	retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  	check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  	retval = IDADenseGetNumResEvals(mem, &nreLS);
  	check_flag(&retval, "IDADenseGetNumResEvals", 1);
  	retval = IDAGetNumGEvals(mem, &nge);
  	check_flag(&retval, "IDAGetNumGEvals", 1);

  	printf("\nFinal Run Statistics: \n\n");
  	printf("Number of steps                    = %ld\n", nst);
  	printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  	printf("Number of Jacobian evaluations     = %ld\n", nje);
  	printf("Number of nonlinear iterations     = %ld\n", nni);
  	printf("Number of error test failures      = %ld\n", netf);
  	printf("Number of nonlinear conv. failures = %ld\n", ncfn);
  	printf("Number of root fn. evaluations     = %ld\n", nge);
  }

/*
 * Check function return value... (Copied from Sundials examples)
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}
