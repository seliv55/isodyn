// clearing vectors
#include <iostream>
#include <vector>

int main ()
{
  std::vector<std::string> myvector;
  myvector.push_back ("asd");
  myvector.push_back ("ert");
  myvector.push_back ("xcv");

  std::cout << "myvector contains:";
  for (unsigned i=0; i<myvector.size(); i++)
    std::cout << ' ' << myvector[i];
  std::cout << '\n';

  myvector.clear();
  myvector.push_back ("abcd");
  myvector.push_back ("efgh");

  std::cout << "myvector contains:";
  for (unsigned i=0; i<myvector.size(); i++)
    std::cout << ' ' << myvector[i];
  std::cout << '\n';

  return 0;
}
