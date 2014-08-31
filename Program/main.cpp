/**************************************************************************\
| File Name: main.cpp                                                      |
| Brief Description: Main-file for Simulating a Solenoid to calculate the  |
| Magnetic field components                                                |
| Date: 12/1/2013                                                          |
| Author: Christopher E. Mertin <cem09v@my.fsu.edu>                        |
| Associated Files: SolenoidBField.h, SolenoidBField.cpp, Variables.h      |
| Detailed Description: This program calculates the off-axis and on-axis   |
| elements of the magnetic field produced from a solenoid. It does this    |
| with the use of the integral form of the Biot-Savart law, with the use of|
| Vector Calculus and numerical integration techniques, instead of         |
| resulting to elliptic integration. This program was made to step through |
| various coordinate values and calculate the magnetic field at each point |
|                                                                          |
| As the calculation can be very time consuming based on the choosen       |
| parameters, the program outputs "progress updates" to standard error so  | 
| that the user knows the current location and progress of the calculations| 
| being performed.                                                         |
|                                                                          |
| The way this program performs the calculations is by treating each loop  |
| of wire in the solenoid as an individual loop of wire, and calculates the|
| magnetic field by summing up the contribution of each loop/coil. This    |
| program was built for the use of simulating multi-layered solenoids,     |
| i.e. solenoids which are wound with more than one layer of wire.         |
\**************************************************************************/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cwchar>

#include "SolenoidBField.h"
#include "Variables.h"

using namespace std;

// Function Declarations ////////////////////////
string convertNum(float number);
string getSysTime();
void calcNorm(struct B_var& objectCalc, struct Coordinate &positions, struct Final &results);
void printHeader(ofstream &paramOutput, struct B_var &Solenoid, struct B_var &Helm1, struct B_var &Helm2, struct Coordinate &positions, struct Final &results, string &runNumber);
void mainCalc(SolenoidBField* B, ofstream &DataOutput, struct B_var &objectCalc, struct Coordinate &positions, struct Final &results);

/******************************************************************\
| Function Name: main()                                            |
| Description: Where all the variables are defined and functions   |
| called in order to calculate the B-field. Determines the way in  |
| which the calculations and integration is done.                  |
\******************************************************************/
int main(int argc, char *argv[])
{
  if(argc != 2)
    {
      cerr << "Illegal number of arguments. Need a single run number." << endl;
      cerr << "Exiting..." << endl << endl;
      return 0;
    }
  string runNumber = argv[1];
  string parametersFilename = "Run" + runNumber + "_Parameters.dat";

  // Solenoid Simulation Constants ////////////////
  B_var Solenoid;
  Solenoid.NIntervals = 10;
  Solenoid.numCoils = 10;//228010; //228010
  Solenoid.length = .63;
  Solenoid.zMin = -(Solenoid.length)/(2.0);
  Solenoid.zMax = (Solenoid.length)/(2.0);
  Solenoid.radiusMin = .14224/2.0;
  Solenoid.radiusMax = .20320/2.0;
  float coilAreaMax = 2 * M_PI * pow(Solenoid.radiusMax,2) + 2 * M_PI * Solenoid.radiusMax * Solenoid.length;
  float coilAreaMin = 2 * M_PI * pow(Solenoid.radiusMin,2) + 2 * M_PI * Solenoid.radiusMin * Solenoid.length;
  float coilArea = coilAreaMax - coilAreaMin;
  float wireArea = coilArea / Solenoid.numCoils;
  Solenoid.wireDiam = sqrt(4 * wireArea / M_PI);
  Solenoid.numCoilAlongZ = Solenoid.length / Solenoid.wireDiam;
  Solenoid.numCoilAlongR = (Solenoid.radiusMax - Solenoid.radiusMin) / Solenoid.wireDiam;
  Solenoid.current = 139.2; //70, 55
  Solenoid.filename = "Solenoid.dat"; //Output file for only Solenoid Data
  Solenoid.xCoil = 0.0; //Initial coil location
  Solenoid.yCoil = 0.0; //Initial coil location
  Solenoid.zCoil = Solenoid.zMax; //Initial coil locatoin
  Solenoid.coilRadius = Solenoid.radiusMin; //Initial radius of coil
  /////////////////////////////////////////////////

  // Helmholtz 1 Simulation Constants /////////////
  B_var Helm1;
  Helm1.NIntervals = 10;
  Helm1.current = 0; //-40;
  Helm1.numCoilAlongZ = 10;
  Helm1.numCoilAlongR = 5;
  Helm1.numCoils = Helm1.numCoilAlongZ * Helm1.numCoilAlongR;
  Helm1.wireDiam = .00145;
  Helm1.length = Helm1.wireDiam * Helm1.numCoilAlongZ;
  Helm1.radiusMin = 0.25;
  Helm1.radiusMax = 0.25 + Helm1.numCoilAlongR * Helm1.wireDiam;
  Helm1.zMin = 1.3 -(Helm1.radiusMin/2.0) - Helm1.length;
  Helm1.zMax = 1.3 -(Helm1.radiusMin/2.0);
  Helm1.filename = "Helmholtz1.dat";
  Helm1.xCoil = 0.0;
  Helm1.yCoil = 0.0;
  Helm1.zCoil = Helm1.zMax;
  Helm1.coilRadius = Helm1.radiusMin;
  /////////////////////////////////////////////////

  // Helmholtz 2 Simulation Constants /////////////
  B_var Helm2;
  Helm2.NIntervals = 10;
  Helm2.current = 0; //50;
  Helm2.numCoilAlongZ = 10;
  Helm2.numCoilAlongR = 5;
  Helm2.numCoils = Helm2.numCoilAlongZ * Helm2.numCoilAlongR;
  Helm2.wireDiam = .00145;
  Helm2.length = Helm2.wireDiam * Helm2.numCoilAlongZ;
  Helm2.radiusMin = 0.25;
  Helm2.radiusMax = 0.25 + Helm2.numCoilAlongR * Helm2.wireDiam;
  Helm2.zMin = 1.3 + (Helm1.radiusMin/2.0);
  Helm2.zMax = 1.3 + (Helm1.radiusMin/2.0) + Helm2.length;
  Helm2.filename = "Helmholtz2.dat";
  Helm2.xCoil = 0.0;
  Helm2.yCoil = 0.0;
  Helm2.zCoil = Helm2.zMax;
  Helm2.coilRadius = Helm2.radiusMin;
  /////////////////////////////////////////////////

  // Positions to Calculate ///////////////////////
  Coordinate positions;
  positions.rPosInit = -0.4;
  positions.rPosFin = 0.4;
  
  positions.zPosInit = -0.6;
  positions.zPosFin = 0.6;
  
  positions.dPos = 0.01;
  positions.xPosInit = positions.rPosInit / sqrt(2.0);
  positions.xPosFin = positions.rPosFin / sqrt(2.0);
  positions.yPosInit = positions.rPosInit / sqrt(2.0);
  positions.yPosFin = positions.rPosFin / sqrt(2.0);

  positions.maxZcalc = (positions.zPosFin - positions.zPosInit)/positions.dPos;
  positions.maxRcalc = (positions.rPosFin - positions.rPosInit)/positions.dPos;
  /////////////////////////////////////////////////

  // Results Parameters and Output ////////////////
  Final results;
  results.numCoordinates = (int)positions.maxZcalc * (1 + (int)positions.maxRcalc);
  results.Bx_tot = 0.0;
  results.By_tot = 0.0;
  results.Bz_tot = 0.0;
  results.Br_tot = 0.0;
  results.B_mag = 0.0;
  results.filename = "Simulation.dat"; 
  ofstream DataOutput; 
  DataOutput.open((results.filename).c_str()); 
  /////////////////////////////////////////////////

  // Array for Calculations ///////////////////////
  results.numRCalc = (int)positions.maxRcalc + 2;
  results.numZCalc = (int)positions.maxZcalc + 2;
  
  results.B_fieldArray_Z = new double*[results.numRCalc]();
  results.B_posArray_Z = new double*[results.numRCalc]();
  results.B_fieldArray_R = new double*[results.numRCalc]();
  results.B_posArray_R = new double*[results.numRCalc]();

  for(int i = 0; i < results.numRCalc; i++)
    {
      results.B_fieldArray_Z[i] = new double[results.numZCalc]();
      results.B_posArray_Z[i] = new double[results.numZCalc]();
      results.B_fieldArray_R[i] = new double[results.numZCalc]();
      results.B_posArray_R[i] = new double[results.numZCalc]();
    }
  /////////////////////////////////////////////////

  // Solenoid Function Calls //////////////////////  
  ofstream paramOutput;
  paramOutput.open(parametersFilename.c_str());

  ofstream DataOutput_S;
  DataOutput_S.open((Solenoid.filename).c_str());
  DataOutput_S << "#   r" << setw(10) << "z" << setw(10) << "Br" << setw(10) << "Bz" << setw(12) << "B_mag" << endl;
  calcNorm(Solenoid, positions, results);
  printHeader(paramOutput, Solenoid, Helm1, Helm2, positions, results, runNumber);
  SolenoidBField* B = new SolenoidBField(Solenoid);
  mainCalc(B, DataOutput_S, Solenoid, positions, results);
  
  DataOutput_S.close();
  results.Bx_tot = results.By_tot = results.Bz_tot = results.Br_tot = results.B_mag = 0;
  /////////////////////////////////////////////////  
  
  // Helmholtz1 Function Calls ////////////////////
  ofstream DataOutput_H1;
  DataOutput_H1.open((Helm1.filename).c_str());
  DataOutput_H1 << "#   r" << setw(10) << "z" << setw(10) << "Br" << setw(10) << "Bz" << setw(12) << "B_mag" << endl;

  SolenoidBField* H1 = new SolenoidBField(Helm1);
  mainCalc(H1, DataOutput_H1, Helm1, positions, results);

  DataOutput_H1.close();
  results.Bx_tot = results.By_tot = results.Bz_tot = results.Br_tot = results.B_mag = 0;
  /////////////////////////////////////////////////  

  // Helmholtz2 Function Calls ////////////////////  
  ofstream DataOutput_H2;
  DataOutput_H2.open((Helm2.filename).c_str());
  DataOutput_H2 << "#   r" << setw(10) << "z" << setw(10) << "Br" << setw(10) << "Bz" << setw(12) << "B_mag" << endl;

  SolenoidBField* H2 = new SolenoidBField(Helm2);
  mainCalc(H2, DataOutput_H2, Helm2, positions, results);

  DataOutput_H2.close();
  /////////////////////////////////////////////////  
  
  float z_hold = 0; //Temporary holder for calculating the final magnitude of the B field
  float r_hold = 0; //Temporary holder for calculating the final magnitude of the B field

  DataOutput << "#  r" << setw(10) << "z" << setw(10) << "Br" << setw(10) << "Bz" << setw(10) << "B_mag" << endl;
  DataOutput << fixed << setprecision(5);
  for(int i = 0; i < results.numRCalc-1; i++)
    {
      for(int j = 0; j < results.numZCalc-1; j++)
	{
	  z_hold = results.B_fieldArray_Z[i][j];
	  r_hold = results.B_fieldArray_R[i][j];
	  DataOutput << results.B_posArray_R[i][j] << setw(10) << results.B_posArray_Z[i][j] << setw(10) << results.B_fieldArray_R[i][j] << setw(10) << results.B_fieldArray_Z[i][j] << setw(10) << sqrt( pow(z_hold,2) + pow(r_hold,2) ) << endl;
	}
      DataOutput << endl;
    }

  DataOutput.close();
  cerr << endl << "Calculations Complete!" << endl;
  /////////////////////////////////////////////////  

  for(int i = 0; i < results.numRCalc; ++i) 
    {
      delete [] results.B_fieldArray_Z[i];
      delete [] results.B_posArray_Z[i];
      delete [] results.B_fieldArray_R[i];
      delete [] results.B_posArray_R[i];
    }
  
  delete [] results.B_fieldArray_Z;
  delete [] results.B_posArray_Z;
  delete [] results.B_fieldArray_R;
  delete [] results.B_posArray_R;

  return 0;
}

/******************************************************************\
| Function Name: convertNum(float number)                          |
| Description: Takes the parameter "number" and converts it to a   |
| string to be used in the filename of the datafile for the        |
| simulation.                                                      |
\******************************************************************/
string convertNum(float number)
{
  stringstream ss; 
  ss << number;
  return ss.str();
}

/******************************************************************\
| Function Name: getSysTime()                                      |
| Description: Pulls the local time from the system and converts   |
| it into a string, and then returns it to where it was called.    |
| This function is used in including the time and date that the    |
| simulation was ran for later referencing. The format of the date |
| and time is as follows:                                          |
|          Simulation started at: 2013-12-17 at 13:37              |
\******************************************************************/
string getSysTime()
{
  int chooser = 0;
  stringstream ss;
  time_t now;
  struct tm *current;
  now = time(0);
  current = localtime(&now);
  ss << "Simulation started at: ";
  for(int i = 0; i < 5; i++)
    {
      switch(chooser)
	{
	case 0:
	  ss << current->tm_year + 1900 << "-"; 
	  break;                                
	case 1:
	  ss << current->tm_mon + 1 << "-"; 
	  break;                            
	case 2:
	  ss << current->tm_mday << " at "; 
	  break;                            
	case 3:
	  ss << current->tm_hour << ":"; 
	  break;                         
	case 4:
	  {
	    if((int)current->tm_min<10)
	      ss << 0 << current->tm_min;
	    else
	      ss << current->tm_min; 
	    break;
	  }                 
	default:
	  ss << "Error" << "-"; 
	  break;                
	}
      chooser++;
    }
  return ss.str();
}

/******************************************************************\
| Function Name: calcNorm(str. B_var, str. Coordinate, str. Final) |
| Description: Calculates the normalization factor for the solenoid|
| which is the stored in the parameters file. The point of this    |
| function is to find the maximum B field of the solenoid, so that |
| value can be used to normalize the largest vector to 1 when      |
| plotting the results.                                            |
\******************************************************************/
void calcNorm(struct B_var &objectCalc, struct Coordinate &positions, struct Final &results)
{
  SolenoidBField* B = new SolenoidBField(objectCalc);

  // Normalization Factor ///////////////////////// 
  cerr << "  * Calculating the Normalization Factor." << endl;
  for(int Rint = 0; Rint < objectCalc.numCoilAlongR; Rint++)   
    {
      objectCalc.coilRadius = objectCalc.radiusMin + Rint * (objectCalc.wireDiam); 
      for(int Zint = 0; Zint < objectCalc.numCoilAlongZ; Zint++) 
	{
	  objectCalc.zCoil = objectCalc.zMax - Zint * (objectCalc.wireDiam); 
	  B->SetCoilDimensions(objectCalc);
	  B->CalculateBField(0.0, 0.0, 0.0);
	  results.Bz_norm += B->GetBz(); 
	}
    }

  delete B;
}

/******************************************************************\
| Function Name: printHeader(ofstream, B_var(3), Coordinate, Final,|
|                            string)                               |
| Description: Creates a file that contains all of the parameters  |
| that were used to run the simulation, as well as the run number  |
| chosen. This is so that when the plots are made, and referenced  |
| at a later time, the run parameters will still be known.         |
\******************************************************************/
void printHeader(ofstream &paramOutput, struct B_var &Solenoid, struct B_var &Helm1, struct B_var &Helm2, struct Coordinate &positions, struct Final &results, string &runNumber)
{
  paramOutput << "# " << setw(50) << "Simulation parameters for Run " << runNumber << endl;
  paramOutput << "# " << setw(55) << getSysTime() << endl; 
  paramOutput << "# " << setw(45) << "Solenoid Values" << endl;
  paramOutput << setprecision(3) << fixed;
  paramOutput << "# Integration Intervals: " << Solenoid.NIntervals << " " << setw(20) << "Number of Coils: " << (int)Solenoid.numCoils << " " << setw(15) << "Length: " << Solenoid.length << endl;
  paramOutput << "# Min z position: " << Solenoid.zMin << setw(25) << "Max z position: " << Solenoid.zMax << endl;
  paramOutput << "# Min Radius: " << Solenoid.radiusMin << setw(15) << "Max Radius: " << Solenoid.radiusMax << setw(13) << "Current: " << Solenoid.current << endl;
  paramOutput << "# Initial Coil Positions: " << endl;
  paramOutput << "# x Coil: " << Solenoid.xCoil << setw(15) << "y Coil: " << Solenoid.yCoil << setw(15) << "z Coil: " << Solenoid.zCoil << endl;
  paramOutput << "# Normalization Factor: " << setw(10) << setprecision(5) << results.Bz_norm << endl;
  paramOutput << setprecision(3) << fixed;
  paramOutput << "# -------------------------------------------------------------------------" << endl;
  paramOutput << "# " << setw(45) << "Helmholtz 1 Values" << endl;
  paramOutput << "# Integration Intervals: " << Helm1.NIntervals << endl;
  paramOutput << "# Number of coils along Z: " << (int)Helm1.numCoilAlongZ << setw(35) << "Number of coils along R: " << (int)Helm1.numCoilAlongR << endl;
  paramOutput << "# Current: " << Helm1.current << setw(25) << "Minimum coil radius: " << Helm1.radiusMin << endl;
  paramOutput << "# Wire Diameter: " << setprecision(5) << Helm1.wireDiam << endl;
  paramOutput << setprecision(3) << fixed;
  paramOutput << "# -------------------------------------------------------------------------" << endl;
  paramOutput << "# " << setw(45) << "Helmholtz 2 Values" << endl;
  paramOutput << "# Integration Intervals: " << Helm2.NIntervals << endl;
  paramOutput << "# Number of coils along Z: " << (int)Helm2.numCoilAlongZ << setw(35) << "Number of coils along R: " << (int)Helm2.numCoilAlongR << endl;
  paramOutput << "# Current: " << Helm2.current << setw(25) << "Minimum coil radius: " << Helm2.radiusMin << endl;
  paramOutput << "# Wire Diameter: " << setprecision(5) << Helm2.wireDiam << endl;
  paramOutput << setprecision(3) << fixed;
  paramOutput << "# -------------------------------------------------------------------------" << endl;
  paramOutput << "# " << setw(45) << "Positions to Calculate" << endl;
  paramOutput << "# Initial x: " << positions.xPosInit << setw(20) << "Initial y: " << positions.yPosInit << setw(20) << "Initial z: " << positions.zPosInit << endl;
  paramOutput << "# Final x: " << positions.xPosFin << setw(21) << "Final y: " << positions.yPosFin << setw(21) << "Final z: " << positions.zPosFin << endl;
  paramOutput << "# Coordinate Precision: " << positions.dPos << endl;
}


/******************************************************************\
| Function Name: mainCalc(SolenoidBField, ofstream, B_var,         |
|                         Coordinate, Final)                       |
| Description: This function performs the main calculation for the |
| simulation by calling the appropriate functions from the         |
| SolenoidBField Class. It stores the data in its own file, and    |
| stores the values in an array so the overall magnetic field can  |
| be found by summing the values later on.                         |
\******************************************************************/
void mainCalc(SolenoidBField* B, ofstream &DataOutput, struct B_var &objectCalc, struct Coordinate &positions, struct Final &results)
{
  // Calculating the Magnetic Field /////////////// 
  float xPos = positions.xPosInit;
  float yPos = positions.yPosInit;
  float zPos = positions.zPosInit;
  float rPos = positions.rPosInit;

  cerr << endl << "& Calculating the magnetic field for " << objectCalc.numCoils << " loops over " << results.numCoordinates << " coordinates." << endl;
  for(int Rcalc = 0; Rcalc <= (int)positions.maxRcalc; Rcalc++) 
    {
      cerr << fixed;
      cerr << endl << setprecision(2) << "& Currently on (z, r): (" << zPos << ", " << rPos << ") - " << Rcalc / positions.maxRcalc * 100 << "% Completed." << endl;
      rPos = positions.rPosInit + Rcalc * positions.dPos; 
      xPos = yPos = rPos/sqrt(2); 

      for(int Zcalc = 0; Zcalc <= (int)positions.maxZcalc; Zcalc++) 
	{
	  results.Bx_tot = 0; 
	  results.By_tot = 0;
	  results.Bz_tot = 0;
	  results.Br_tot = 0;
	  results.B_mag = 0;
	  zPos = positions.zPosInit + Zcalc * positions.dPos; 
	  for(int Rint = 0; Rint < objectCalc.numCoilAlongR; Rint++) 
	    {
	      objectCalc.coilRadius = objectCalc.radiusMin + Rint * objectCalc.wireDiam; 
	      for(int Zint = 0; Zint < objectCalc.numCoilAlongZ; Zint++) 
		{
		  objectCalc.zCoil = objectCalc.zMax - Zint * objectCalc.wireDiam; 
		  B->SetCoilDimensions(objectCalc);
		  B->CalculateBField(xPos, yPos, zPos);
		  results.Bx_tot += B->GetBx(); 
		  results.By_tot += B->GetBy(); 
		  results.Bz_tot += B->GetBz(); 
		}
	    }
	  results.Br_tot = sqrt(pow(results.Bx_tot,2) + pow(results.By_tot,2)); 
	  if(objectCalc.current < 0) //Changes direction based on current
	    {
	      results.Br_tot = -results.Br_tot;
	    }
	  if((rPos <= 0 && zPos >= 0) || (rPos >= 0 && zPos <= 0)) //Makes B field symmetric 
	    {
	      results.Br_tot = -results.Br_tot;
	    }
	  results.B_mag = sqrt(pow(results.Br_tot,2) + pow(results.Bz_tot,2)); 
	  DataOutput << fixed;
	  DataOutput << setprecision(5) << rPos << setw(10) << zPos << setw(10) << results.Br_tot << setw(10) << results.Bz_tot << setw(10) << results.B_mag << endl; 
	  results.B_fieldArray_Z[Rcalc][Zcalc] += results.Bz_tot;
	  results.B_posArray_Z[Rcalc][Zcalc] = zPos;
	  results.B_fieldArray_R[Rcalc][Zcalc] += results.Br_tot;
	  results.B_posArray_R[Rcalc][Zcalc] = rPos;

	}
      DataOutput << endl;
    }
}
