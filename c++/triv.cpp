#include <iostream>
#include <string>


int main (int argc, char** argv){
	std::string s1, s2;
	s1 = "boots";
	s2 = "swag";

	std::cout<<s1 <<" "<<s2;

	s1 = s2;
	s2 = "gimli";
	
	
	std::cout<<s1 <<" "<<s2;
	return 0;
}
