#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;
int count(int ini=1);

int count(int ini){
    int ifin;
    char fn[5];
    for(int i=ini;;i++) { sprintf(fn,"%i",i);
	   ifstream checkfi(fn);
	   if(!checkfi.good()) { ifin=i-1; break;}
	   checkfi.close();
   }
   return ifin;}
   
void fimove(int inew, int init, int ifin){
   char fn[15]; 
  for (int i=init;i<ifin;i++) {sprintf(fn,"cp %i ../files/%i",i,(i+inew-init));
		system(fn);}
}
    
int main( int argc, char *argv[] ){
   int init, ifin, inew;
   cout<<"Old first file="; cin>>init;
   cout<<"New first file="; cin>>inew;
	 ifin=count(init);
	 
	 fimove(inew,init,ifin);
    cout <<"init="<<init<<"; ifin=" <<ifin<<"; inew=" <<inew<<endl;
return 0;
}
