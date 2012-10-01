h = 0.00625;

Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};
Point(3) = {1,1,0,h};
Point(4) = {0,1,0,h};
Point(5) = {0,0,1,h};
Point(6) = {1,0,1,h};
Point(7) = {1,1,1,h};
Point(8) = {0,1,1,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {1,5};
Line(6) = {2,6};
Line(7) = {3,7};
Line(8) = {4,8};
Line(9) = {5,6};
Line(10) = {6,7};
Line(11) = {7,8};
Line(12) = {8,5};

Line Loop(1) = {-1,-4,-3,-2};
Line Loop(2) = {1,6,-9,-5};
Line Loop(3) = {2,7,-10,-6};
Line Loop(4) = {3,8,-11,-7};
Line Loop(5) = {4,5,-12,-8};
Line Loop(6) = {9,10,11,12};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop (1) = {1,2,3,4,5,6};

Volume (1) = {1};
