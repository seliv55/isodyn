#include <iostream>
using namespace std;

int main ( int argc, char *argv[] ) {
//const char* b=a.c_str();
//x += atoi(y);  // Use the integer equivalent of char / string  y;
 int b[]={1,2,5,18}, *ii[2]; ii[0]=&b[1]; ii[1]=&b[1];
 for (int i=0;i<4;i++) b[i] += i;
 int i=((31>>1)&((32>>1>>1)-1));
cout << ii[0][2]<<"  "<<ii[1][2] <<" i="<<i<< endl;
return 0;}
