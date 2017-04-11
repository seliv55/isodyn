#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

int main() {
    using namespace std;
    string sentence = "And/I/feel/fine, but...";
    istringstream iss(sentence);
    stringstream jss(sentence);
    string sss;
    while(getline(jss,sss,'/')) cout<<sss<<"\n";
    int i=sentence.find_last_of('/');  cout<<i<<"\n";
    string a=sentence.substr(0,i);  cout<<a<<"\n";
    string b=sentence.substr(i+1);  cout<<b<<"\n";
}
