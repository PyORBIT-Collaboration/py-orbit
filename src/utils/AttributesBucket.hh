//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   AttributesBucket.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    12/14/2005
//
// DESCRIPTION
//    A class for collection of attributes' . This is a class where user keeps
//    different data in the form key-value.
//
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#ifndef ATTRIBUTES_BUCKET_H
#define ATTRIBUTES_BUCKET_H

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include <string>
#include <map>
#include <vector>

#include <iostream>
#include <fstream>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    AttributesBucket
//
///////////////////////////////////////////////////////////////////////////

class AttributesBucket
{
  //--------------------------------------
  //the public methods of the AttributesBucket class
  //--------------------------------------
public:

  AttributesBucket();
  ~AttributesBucket();

  int intVal(const std::string attName);
  double doubleVal(const std::string attName);
	
  int intVal(const std::string attName, int val);
  double doubleVal(const std::string attName, double val);

  //returns 0 if there is no attribute associated with this name
  int hasIntAttribute(const std::string attName);
  int hasDoubleAttribute(const std::string attName);

  void getIntAttributeNames(std::vector<std::string>& names);
  void getDoubleAttributeNames(std::vector<std::string>& names);

	//add all atributes from one bucket to another
	void add(AttributesBucket* bckt);

	//removes all values
	void clear();

  //-----------------------------------
  //  DATA MEMBERS
  //-----------------------------------

private:

  std::map<std::string,int> intAttrMap;
  std::map<std::string,double> doubleAttrMap;

};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
