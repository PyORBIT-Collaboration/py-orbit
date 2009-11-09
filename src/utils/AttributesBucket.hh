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
//    A class for collection of int and double attributes in the dictionary. 
///////////////////////////////////////////////////////////////////////////

#ifndef ATTRIBUTES_BUCKET_H
#define ATTRIBUTES_BUCKET_H

#include <string>
#include <map>
#include <vector>

namespace OrbitUtils{ 
	
/** 
  A class is a collection of int and double key(string)-value pairs. 
  Users have to keep in mind that operations with dictionaries are slow. 
*/


	class AttributesBucket
	{
		//--------------------------------------
		//the public methods of the AttributesBucket class
		//--------------------------------------
	public:
		
		/** Constructor for an empty dictinary. */
		AttributesBucket();
		
		/** Destructor. */
		~AttributesBucket();
		
		/** Returns the int value from the dictionary for the key string. */
		int intVal(const std::string attName);
		
		/** Returns the double value from the dictionary for the key string. */
		double doubleVal(const std::string attName);
		
		/** Sets the int value to the dictionary for the key string. */
		int intVal(const std::string attName, int val);
		
		/** Sets the double value to the dictionary for the key string. */
		double doubleVal(const std::string attName, double val);
		
		/** Returns 0 if there is no attribute associated with this name */
		int hasIntAttribute(const std::string attName);
		
		/** Returns 0 if there is no attribute associated with this name */
		int hasDoubleAttribute(const std::string attName);
		
		/** Returns a vector with integer attributes keys. */
		void getIntAttributeNames(std::vector<std::string>& names);
		
		/** Returns a vector with doule attributes keys. */
		void getDoubleAttributeNames(std::vector<std::string>& names);
		
		/** Adds all atributes from one bucket to another. */
		void addTo(AttributesBucket* bckt);
		
		/** Removes all values from the dictionary. */
		void clear();
		
		//-----------------------------------
		//  DATA MEMBERS
		//-----------------------------------
		
	private:
		
		std::map<std::string,int> intAttrMap;
		std::map<std::string,double> doubleAttrMap;
		
	};

};

#endif
