/**************************************************************************\
| File Name: SolenoidBField.h                                              |
| Brief Description: Header file for Simulating a Solenoid to calculate the|
| Magnetic field components                                                |
| Date: 12/1/2013                                                          |
| Author: Christopher E. Mertin <cem09v@my.fsu.edu>                        |
| Associated Files: main.cpp & SolenoidBField.cpp                          |
| Detailed Description: For use with the main.cpp and SolenoidBField.cpp   |
| files that were made along side this file, with the purpose of           |
| calculating the radial and axial components of a magnetic field produced |
| by a solenoid. This SolenoidBField class holds the funciton declarations |
| which are defined in the SolenoidBField.cpp file, which is the main      |
| source of the magnetic field calculations.                               |
\**************************************************************************/

#ifndef SOLENOID_B_FIELD_H
#define SOLENOID_B_FIELD_H

#include <iostream>
#include "Variables.h"

using namespace std;

class SolenoidBField
{
 public:
  SolenoidBField(struct B_var &objectCalc);
  void CalculateBField(float xPos, float yPos, float zPos); 
  void SetCoilDimensions(struct B_var &objectCalc);
  //float xComp(float x, float y, float z); 
  //float yComp(float x, float y, float z); 
  //float zComp(float x, float y, float z); 
  int GetNIntervals(); 
  float GetRadius(); 
  float GetCurrent(); 
  float GetBx(); 
  float GetBy(); 
  float GetBz(); 
  float GetxCoil(); 
  float GetyCoil(); 
  float GetzCoil();
  float Bx_fn(float *xyzPos, float t);
  float By_fn(float *xyzPos, float t);
  float Bz_fn(float *xyzPos, float t);
  float Integrate(float (SolenoidBField::*)(float * array, float t), float * xyzPos);

 private:
  // Magnetic Variables
  int NIntervals;
  float numCoils;
  float length;
  float zMin;
  float zMax;
  float radiusMin;
  float radiusMax;
  float wireDiam;
  float numCoilAlongZ;
  float numCoilAlongR;
  float current;
  string filename;
  ofstream DataOutput;

  // Changing Magnetic Variables
  float radius;
  float xCoil, yCoil, zCoil;
  float Bx, By, Bz;
  float Bx_tot, By_tot, Bz_tot;
  float Bz_norm;
  float B_mag;

  // Position Variables
  float rPosInit;
  float rPosFin;
  float zPosInit;
  float zPosFin;
  float dPos;
  float xPosInit;
  float xPosFin;
  float yPosInit;
  float yPosFin;
  float maxZcalc;
  float maxRcalc;
};

#endif
