#include <iostream>
#include <limits>

#include "Rivet/Tools/Cmp.hh"

using namespace std;

int main() {
  using namespace Rivet;

  CmpState cs = CmpState::UNDEF;

  cs = cmp(0.5, 0.6);
  cout << "cmp(0.5, 0.6) = " << cs << '\n';
  assert(cs == CmpState::NEQ);

  cs = cmp(0.5, 0.5);
  cout << "cmp(0.5, 0.5) = " << cs << '\n';
  assert(cs == CmpState::EQ);

  cs = cmp(0.6, 0.5);
  cout << "cmp(0.6, 0.5) = " << cs << '\n';
  assert(cs == CmpState::NEQ);

  cs = cmp(1.,1.) || cmp(0.6, 0.5);
  cout << "cmp(1.,1.) || cmp(0.6, 0.5) = " << cs << '\n';
  assert(cs == CmpState::NEQ);

  return EXIT_SUCCESS;
}
