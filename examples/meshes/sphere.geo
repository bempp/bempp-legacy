h = 0.1;

Point(1) = {0,0,0,h};
Point(2) = {1,0,0,h};
Point(3) = {0,1,0,h};
Point(4) = {0,0,1,h};
Point(5) = {-1,0,0,h};
Point(6) = {0,-1,0,h};
Point(7) = {0,0,-1,h};

Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,6};
Circle(4) = {6,1,2};
Circle(5) = {2,1,7};
Circle(6) = {7,1,5};
Circle(7) = {5,1,4};
Circle(8) = {4,1,2};
Circle(9) = {6,1,7};
Circle(10) = {7,1,3};
Circle(11) = {3,1,4};
Circle(12) = {4,1,6};

Line Loop(1) = {1,11,8};
Line Loop(2) = {2,7,-11};
Line Loop(3) = {3,-12,-7};
Line Loop(4) = {4,-8,12};
Line Loop(5) = {5,10,-1};
Line Loop(6) = {-2,-10,6};
Line Loop(7) = {-3,-6,-9};
Line Loop(8) = {-4,9,-5}; 

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};
Ruled Surface(7) = {7};
Ruled Surface(8) = {8};

Surface Loop (1) = {1,2,3,4,5,6,7,8};

Volume (1) = {1};

Point(11) = {0,3,0,h};
Point(12) = {1,3,0,h};
Point(13) = {0,4,0,h};
Point(14) = {0,3,1,h};
Point(15) = {-1,3,0,h};
Point(16) = {0,-2,0,h};
Point(17) = {0,3,-1,h};

Circle(21) = {12,11,13};
Circle(22) = {13,11,15};
Circle(23) = {15,11,16};
Circle(24) = {16,11,12};
Circle(25) = {12,11,17};
Circle(26) = {17,11,15};
Circle(27) = {15,11,14};
Circle(28) = {14,11,12};
Circle(29) = {16,11,17};
Circle(210) = {17,11,13};
Circle(211) = {13,11,14};
Circle(212) = {14,11,16};

Line Loop(31) = {21,211,28};
Line Loop(32) = {22,27,-211};
Line Loop(33) = {23,-212,-27};
Line Loop(34) = {24,-28,212};
Line Loop(35) = {25,210,-21};
Line Loop(36) = {-22,-210,26};
Line Loop(37) = {-23,-26,-29};
Line Loop(38) = {-24,29,-25}; 

Ruled Surface(41) = {31};
Ruled Surface(42) = {32};
Ruled Surface(43) = {33};
Ruled Surface(44) = {34};
Ruled Surface(45) = {35};
Ruled Surface(46) = {36};
Ruled Surface(47) = {37};
Ruled Surface(48) = {38};

Surface Loop (2) = {41,42,43,44,45,46,47,48};

Volume (2) = {2};