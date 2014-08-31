#ifndef VARIABLES_H
#define VARIABLES_H

#include <iostream>
#include <fstream>

using namespace std;

struct B_var
{
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
  float xCoil;
  float yCoil;
  float zCoil;
  float coilRadius;
  string filename;
  ofstream DataOutput;
};

struct Coordinate
{
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

struct Final
{
  int numCoordinates;
  int numRCalc;
  int numZCalc;
  float normConst;
  float Bx_tot;
  float By_tot;
  float Bz_tot;
  float Br_tot;
  float Bz_norm;
  float B_mag;
  string filename;
  double **B_fieldArray_Z;
  double **B_fieldArray_R;
  double **B_posArray_Z;
  double **B_posArray_R;
};

#endif
