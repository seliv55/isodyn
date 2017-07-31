#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;
    
int main( int argc, char *argv[] ){
   int iflast=2;
   string fidir;
   cout<<"directory="; cin>>fidir;
   cout<<"How many files to keep="; cin>>iflast;
for(int i=iflast+1;;i++) {std::stringstream fn;
       fn<<fidir<<i;   std::ifstream checkfi(fn.str().c_str()); checkfi.close(); 
	 if(checkfi.good()){
	 string com="rm "+ fn.str();  int sys=system(com.c_str());
	 } else return 0;}
	
return 0;
}
