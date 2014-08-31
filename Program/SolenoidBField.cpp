/**************************************************************************\
| File Name: SolenoidBField.cpp                                            |
| Brief Description: File which holds the function definitions in          |
| calculating the magnetic field of a solenoid.                            |
| Date: 12/1/2013                                                          |
| Author: Christopher E. Mertin <cem09v@my.fsu.edu>                        |
| Associated Files: SolenoidBField.h & main.cpp                            |
| Detailed Description: For use with the main.cpp and SolenoidBField.h     |
| files that were made along side this file, with the purspose of          |
| calculating th radial and axial components produced by a multi-layered   |
| solenoid. This file defines all of the function calls for the            |
| SolenoidBField class, which is defined in the SolenoidBField.h file.     |
\**************************************************************************/

#include <iostream>
#include <cmath>

#include "SolenoidBField.h"
#include "Variables.h"

using namespace std;


/******************************************************************\
| Function Name: SolenoidBField(int, float, float, float, float)   |
| Description: Constructor for the SolenoidBField class, it sets   |
| all the passed in parameters as object values which are defined  |
| \listed in the header file.                                      |
\******************************************************************/
SolenoidBField::SolenoidBField(struct B_var &objectCalc)
{
  this->NIntervals = objectCalc.NIntervals;
  this->numCoils = objectCalc.numCoils;
  this->length = objectCalc.length;
  this->zMin = objectCalc.zMin;
  this->zMax = objectCalc.zMax;
  this->radiusMin = objectCalc.radiusMin;
  this->radiusMax = objectCalc.radiusMax;
  this->wireDiam = objectCalc.wireDiam;
  this->numCoilAlongZ = objectCalc.numCoilAlongZ;
  this->numCoilAlongR = objectCalc.numCoilAlongR;
  this->current = objectCalc.current;
  this->filename = objectCalc.filename;
}

/******************************************************************\
| Function Name: CalculateBField(float, float, float)              |
| Description: Calls all the appropriate functions for calculating |
| the x, y, and z components of the B-field. Also separates the    |
| constants, making it easier to folllow the math.                 |
\******************************************************************/
void SolenoidBField::CalculateBField(float xPos, float yPos, float zPos)
{
  float permConst = 4 * M_PI * pow(10, -7); 
  float frontConstant = (permConst * GetCurrent() * GetRadius())/(4 * M_PI); 
  float tempBx, tempBy, tempBz;
  float *xyzPos = new float[3];
  xyzPos[0] = xPos;
  xyzPos[1] = yPos;
  xyzPos[2] = zPos;
  tempBx = frontConstant * (zPos - GetzCoil()) * Integrate(&SolenoidBField::Bx_fn, xyzPos);
  tempBy = frontConstant * (zPos - GetzCoil()) * Integrate(&SolenoidBField::By_fn, xyzPos);  
  tempBz = frontConstant * Integrate(&SolenoidBField::Bz_fn, xyzPos);

  this->Bx = tempBx;
  this->By = tempBy;
  this->Bz = tempBz;  

  delete[] xyzPos;
}

/******************************************************************\
| Function Name: Bx_fn(float *, float)                             |
| Description: Function for returning the value of the x component |
| to the Integrator function, which is integrating over these      |
| values.                                                          |
\******************************************************************/
float SolenoidBField::Bx_fn(float *xyzPos, float t)
{
  float tempZ = xyzPos[2] - GetzCoil();
  float tempR = GetRadius();
  float Denominator = pow(xyzPos[0],2) + pow(xyzPos[1],2) + pow(tempR,2) - 
    2 * tempR * (xyzPos[0] * cos(t) + xyzPos[1] * sin(t)) + pow(tempZ,2);
  Denominator = pow(Denominator,1.5);
  return (cos(t) / Denominator);
}

/******************************************************************\
| Function Name: By_fn(float *, float)                             |
| Description: Function for returning the value of the y component |
| to the Integrator function, which is integrating over these      |
| values.                                                          |
\******************************************************************/
float SolenoidBField::By_fn(float *xyzPos, float t)
{
  float tempZ = xyzPos[2] - GetzCoil();
  float tempR = GetRadius();
  float Denominator = pow(xyzPos[0],2) + pow(xyzPos[1],2) + pow(tempR,2) - 2 * tempR * (xyzPos[0] * cos(t) + xyzPos[1] * sin(t)) + pow(tempZ,2);
  Denominator = pow(Denominator,1.5);
  return (sin(t) / Denominator);
}

/******************************************************************\
| Function Name: Bz_fn(float *, float)                             |
| Description: Function for returning the value of the z component |
| to the Integrator function, which is integrating over these      |
| values.                                                          |
\******************************************************************/
float SolenoidBField::Bz_fn(float *xyzPos, float t)
{
  float tempZ = xyzPos[2] - GetzCoil();
  float tempR = GetRadius();
  float Numerator = tempR - xyzPos[1] * sin(t) - xyzPos[0] * cos(t);
  float Denominator = pow(xyzPos[0],2) + pow(xyzPos[1],2) + pow(tempR,2) - 2 * tempR * (xyzPos[0] * cos(t) + xyzPos[1] * sin(t)) + pow(tempZ,2);
  Denominator = pow(Denominator,1.5);
  return (Numerator / Denominator);
}

/******************************************************************\
| Function Name: Integrate(float (SolenoidBField::*), float *)     |
| Description: Integrates over the functions Bx_fn, By_fn, and     |
| Bz_fn by using the adaptive quadrature method. This method was   |
| taken from "Numerical Analysis 8th Edition" by Burden.           |
\******************************************************************/
float SolenoidBField::Integrate(float (SolenoidBField::*Bxyz)(float *array, float t), float *xyzPos)
{
  int tempN = GetNIntervals();
  float low = 0.0;
  float high = 2.0 * M_PI;
  float TOL_max = pow(10,-5);
  int i = 1;
  float APP = 0.0;
  float FD = 0.0;
  float FE = 0.0;
  float S1 = 0.0;
  float S2 = 0.0;
  int * L = new int[tempN + 1];
  float * v = new float[10];
  float * a = new float[tempN + 1];
  float * h = new float[tempN + 1];
  float * S = new float[tempN + 1];
  float * FA = new float[tempN + 1];
  float * FB = new float[tempN + 1];
  float * FC = new float[tempN + 1];
  float * TOL = new float[tempN + 1];

  TOL[i] = 10.0 * TOL_max;
  a[i] = low;
  h[i] = (high - low)/2.0;
  FA[i] = (this->*Bxyz)(xyzPos, a[i]);
  FC[i] = (this->*Bxyz)(xyzPos, a[i] + h[i]);
  FB[i] = (this->*Bxyz)(xyzPos, high);
  S[i] = h[i] * (FA[i] + 4 * FC[i] + FB[i])/3.0; // Simpson's Rule Approximation over entire interval
  L[i] = 1;

  while(i > 0)
    {
      FD = (this->*Bxyz)(xyzPos, a[i] + h[i]/2);
      FE = (this->*Bxyz)(xyzPos, a[i] + 3 * h[i]/2);
      S1 = h[i] * (FA[i] + 4 * FD + FC[i]) / 6; // Approximations from Simpson's for halves of subintervals
      S2 = h[i] * (FC[i] + 4 * FE + FB[i]) / 6;
      v[1] = a[i];
      v[2] = FA[i];
      v[3] = FC[i];
      v[4] = FB[i];
      v[5] = h[i];
      v[6] = TOL[i];
      v[7] = S[i];
      v[8] = L[i];

      i--;
      if( abs(S1 + S2 - v[7]) < v[6] )
	{
	  APP += S1 + S2;
	}
      else
	{
	  if((int)v[8] >= tempN)
	    {
	      //if(xyzPos[0] >= 0.00005 || xyzPos[1] >= 0.00005)
	      //cerr << "Interval Level Exceeded!" << endl;
	      break;
	    }
	  // Data for right half subinterval
	  i++; 
	  a[i] = v[1] + v[5];
	  FA[i] = v[3];
	  FC[i] = FE;
	  FB[i] = v[4];
	  h[i] = v[5]/2.0;
	  TOL[i] = v[6]/2.0;
	  S[i] = S2;
	  L[i] = v[8] + 1;
	  // Data for left half subinterval
	  i++;
	  a[i] = v[1];
	  FA[i] = v[2];
	  FC[i] = FD;
	  FB[i] = v[3];
	  h[i] = h[i-1];
	  TOL[i] = TOL[i-1];
	  S[i] = S1;
	  L[i] = L[i-1];
	}
    }
  delete[] v;
  delete[] a;
  delete[] h;
  delete[] S;
  delete[] L;
  delete[] FA;
  delete[] FB;
  delete[] FC;
  delete[] TOL;

  return APP;
}

/******************************************************************\
| Function Name: SetCoilDimensions(float, float, float, float)     |
| Description: Sets the dimensions of the current loop of wire/coil|
| to be used in the next B-field calculation.                      |
\******************************************************************/
void SolenoidBField::SetCoilDimensions(struct B_var &objectCalc)
{
  this->xCoil = xCoil;
  this->yCoil = yCoil;
  this->zCoil = zCoil;
  this->radius = objectCalc.coilRadius;
}

int SolenoidBField::GetNIntervals()
{
  return NIntervals;
}

float SolenoidBField::GetRadius()
{
  return radius;
}

float SolenoidBField::GetCurrent()
{
  return current;
}

float SolenoidBField::GetBx()
{
  return Bx;
}

float SolenoidBField::GetBy()
{
  return By;
}

float SolenoidBField::GetBz()
{
  return Bz;
}

float SolenoidBField::GetxCoil()
{
  return xCoil;
}

float SolenoidBField::GetyCoil()
{
  return yCoil;
}

float SolenoidBField::GetzCoil()
{
  return zCoil;
}

