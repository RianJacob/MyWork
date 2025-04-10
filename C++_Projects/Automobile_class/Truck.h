#ifndef Truck_H
#define Truck_H
#include "Automobile.h"
#include <string>
using namespace std;

// The Truck class represents a truck.
class Truck : public Automobile
{
 private:
  string driveType;

 public:
  // Default constructor
 Truck() : Automobile()
    {  driveType = ""; }

  // Constructor #2
 Truck(string truckMake, int truckModel, int truckMileage, double truckPrice, string truckDriveType) :
  Automobile(truckMake, truckModel, truckMileage, truckPrice)
    { driveType = truckDriveType; }

  // Accessor for riveType attribute
  string getDriveType()
  { return driveType; }
};

#endif
  
