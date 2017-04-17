// string::substr
#include <iostream>
#include <string>

int main ()
{
  char *st="We think in generalities, but we live in details.";
   std::string str(st);
                                          // (quoting Alfred N. Whitehead)

  std::string str2 = str.substr (3,5);     // "think"

  std::size_t pos = str.find("but");      // position of "live" in str

  std::string str3 = str.substr (pos);     // get from "live" to the end

  std::cout << str2 << ' ' << str3 << '\n';

  return 0;
}
