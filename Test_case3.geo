

lc = 0.02; //0.005;
Point(1) = {0,-0.15,0,lc}; //{0,-0.01,0,lc};
Point(2) = {1,-0.15,0,lc}; //{1,-0.01,0,lc};
Point(3) = {1,0.15,0,lc}; //{1,0.01,0,lc};
Point(4) = {0,0.15,0,lc}; //{0,0.01,0,lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
ll = newll; Curve Loop(newll) = {1,2,3,4};
Plane Surface(1) = {ll};
Physical Curve("bottom", 1) = {1};
Physical Curve("right", 2) = {2};
Physical Curve("top", 3) = {3};
Physical Curve("left", 4) = {4};
Physical Curve("Bnd_D", 50) = {1,2,3,4};
Physical Surface("D", 100) = {1};



