#ifndef __STRING_UTILS_H_
#define __STRING_UTILS_H_

#include <string>
#include <vector>

using namespace std;

//=======================================
//Usage
// int nTokens = StringUtils::Tokenize("....",vector<string>& tokens, ", :<>");
//or for space-delimiter
//int nTokens = StringUtils::Tokenize("....",vector<string>& tokens);
//=======================================


namespace OrbitUtils{
	namespace StringUtils
	{
		int Tokenize(const string& str,vector<string>& tokens,  const string& delimiters = " ");
	};
};
#endif

