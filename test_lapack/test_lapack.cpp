/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <complex>
#include <iostream>

extern "C" {
  void ilaver_(int* vers_major, int* vers_minor, int* vers_patch);
}

int main()
{
  int vers_major;
  int vers_minor;
  int vers_patch;
  ilaver_(&vers_major, &vers_minor, &vers_patch);
  std::cout << "LAPACK Version ";
  std::cout << vers_major << "." << vers_minor << "." << vers_patch
	    << " detected" << std::endl;

  return 0;
}
