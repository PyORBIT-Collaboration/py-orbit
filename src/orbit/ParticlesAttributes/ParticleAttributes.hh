//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ParticleAttributes.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    07/14/2005
//
// DESCRIPTION
//    An abstract class for particle attributes. This is a base class for the
//    others classes that keep different attributes, e.g. spin orientation, tunes,
//    quantum amplitudes of exited states etc.
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#ifndef PARTICLE_ATTRIBUTES_H
#define PARTICLE_ATTRIBUTES_H

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include <string>
#include <map>

///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    ParticleAttributes
//
///////////////////////////////////////////////////////////////////////////

//forward definition
class  Bunch;

class ParticleAttributes
{
public:
  //--------------------------------------
  //the public methods of the ParticleAttributes class
  //--------------------------------------

  ParticleAttributes(Bunch* bunch, int size_in);
  virtual ~ParticleAttributes();

  //returns the name of the particle attributes bucket
  const std::string& name();

  //returns the attributes description for this particular attributes
  //order and dimensions
  const std::string& attrDescription();

  //returns the pointer to bunch
  const Bunch* bunch();

  //returns the attribute value with particular index for particular particle
  double& attValue(int particle_index, int att_index);

  //returns the attribute array for particular particle
  double* attArr(int particle_index);

  //returns the size of the attrubute bucket
  virtual int getAttSize();

  //returns the number of particles
  int getBunchSize();

  //returns the flag of the particle (dead or alive)
  int flag(int particle_index);

public:
	
	//dictionary with parameters. It is used during the dump or read procedures.
	std::map<std::string,double> parameterDict;
	
protected:

  //initializes the attribute data
  virtual void init(int particle_index);

  //returns the attributes shift from
  //the beginning the combined attribute array in the Bunch
  int getAttrShift();

  //sets the attributes shift from
  //the beginning the combined attribute array in the Bunch
  void setAttrShift(int attr_ind_shift);

  //class bunch can access to the protected methods
  //of the ParticleAttributes class
  friend class Bunch;

  //-----------------------------------
  //  DATA MEMBERS
  //-----------------------------------

  std::string cl_name_;

  std::string attrDescr;

  Bunch* bunch_;

  int attr_ind_shift_;
	
	//size of the array for particles attribute
	int size;

private:


};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
