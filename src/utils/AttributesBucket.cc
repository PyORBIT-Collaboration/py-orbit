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

#include "AttributesBucket.hh"

using namespace OrbitUtils;

AttributesBucket::AttributesBucket()
{
}

AttributesBucket::~AttributesBucket()
{
}
int AttributesBucket::intVal(const std::string attName){
	if(hasIntAttribute(attName) == 0){
		return 0;
	}
	return intAttrMap[attName];
}

int AttributesBucket::intVal(const std::string attName, int val){
	intAttrMap[attName] = val;
	return val;
}

double AttributesBucket::doubleVal(const std::string attName){
	if(hasDoubleAttribute(attName) == 0){
		return 0.;
	}
	return doubleAttrMap[attName];
}

double AttributesBucket::doubleVal(const std::string attName, double val){
	doubleAttrMap[attName] = val;
	return val;
}


int AttributesBucket::hasIntAttribute(const std::string attName){
  return intAttrMap.count(attName);
}

int AttributesBucket::hasDoubleAttribute(const std::string attName){
  return doubleAttrMap.count(attName);
}

void AttributesBucket::getIntAttributeNames(std::vector<std::string>& names){
  names.clear();
  std::map<std::string,int>::iterator pos;
  for (pos = intAttrMap.begin(); pos != intAttrMap.end(); ++pos) {
    std::string name = pos->first;
    names.push_back(name);
  }
}

void AttributesBucket::getDoubleAttributeNames(std::vector<std::string>& names){
  names.clear();
  std::map<std::string,double>::iterator pos;
  for (pos = doubleAttrMap.begin(); pos != doubleAttrMap.end(); ++pos) {
    std::string name = pos->first;
    names.push_back(name);
  }
}

void AttributesBucket::addTo(AttributesBucket* bckt){
	std::vector<std::string> names;
	bckt->getIntAttributeNames(names);
	for(int i = 0, n = names.size(); i < n; i++){
		bckt->intVal(names[i],intVal(names[i]));
	}
	
	names.clear();
	bckt->getDoubleAttributeNames(names);
	for(int i = 0, n = names.size(); i < n; i++){
		bckt->doubleVal(names[i],doubleVal(names[i]));
	}
}

void AttributesBucket::clear(){
	intAttrMap.clear();
	doubleAttrMap.clear();
}


