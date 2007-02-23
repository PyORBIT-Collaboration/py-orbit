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
#include "AttributesBucket.hh"

///////////////////////////////////////////////////////////////////////////
//   Constructor and Desctructor
///////////////////////////////////////////////////////////////////////////

AttributesBucket::AttributesBucket()
{

}

AttributesBucket::~AttributesBucket()
{

}

int& AttributesBucket::intVal(const std::string attName){
  return intAttrMap[attName];
}

double& AttributesBucket::doubleVal(const std::string attName){
   return doubleAttrMap[attName];
}

//returns true if it has int attribute with this key
int AttributesBucket::hasIntAttribute(const std::string attName){
  return intAttrMap.count(attName);
}

//returns true if it has double attribute with this key
int AttributesBucket::hasDoubleAttribute(const std::string attName){
  return doubleAttrMap.count(attName);
}

//returns a vector with integer attributes keys
void AttributesBucket::getIntAttributeNames(std::vector<std::string>& names){
  names.clear();
  std::map<std::string,int>::iterator pos;
  for (pos = intAttrMap.begin(); pos != intAttrMap.end(); ++pos) {
    std::string name = pos->first;
    names.push_back(name);
  }
}

//returns a vector with doule attributes keys
void AttributesBucket::getDoubleAttributeNames(std::vector<std::string>& names){
  names.clear();
  std::map<std::string,double>::iterator pos;
  for (pos = doubleAttrMap.begin(); pos != doubleAttrMap.end(); ++pos) {
    std::string name = pos->first;
    names.push_back(name);
  }
}

//add all atributes from one bucket to another
void AttributesBucket::add(AttributesBucket* bckt){
	std::vector<std::string> names;
	 bckt->getIntAttributeNames(names);
	 for(int i = 0, n = names.size(); i < n; i++){
		 intVal(names[i]) = bckt->intVal(names[i]);
	 }

	 names.clear();
	 bckt->getDoubleAttributeNames(names);
	 for(int i = 0, n = names.size(); i < n; i++){
		 doubleVal(names[i]) = bckt->doubleVal(names[i]);
	 }
}

//removes all values
void AttributesBucket::clear(){
	intAttrMap.clear();
	doubleAttrMap.clear();
}


