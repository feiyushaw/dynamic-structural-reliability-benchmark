// Cantilever Steel Beam Geometry for Dynamic Structural Reliability Validation
// Dimensions: 300 mm (Length) x 40 mm (Width/Height) x 10 mm (Thickness)
// All units in meters

L = 0.300; // Length (X)
H = 0.040; // Height (Y)
B = 0.010; // Thickness (Z)

lc = 0.005; // Target mesh size (5 mm)

// Define Points
Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, H, 0, lc};
Point(4) = {0, H, 0, lc};

Point(5) = {0, 0, B, lc};
Point(6) = {L, 0, B, lc};
Point(7) = {L, H, B, lc};
Point(8) = {0, H, B, lc};

// Define Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Define Surfaces
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1}; // Back Face (Z = 0)

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2}; // Front Face (Z = B)

Curve Loop(3) = {4, 9, -8, -12};
Plane Surface(3) = {3}; // Fixed End Face (X = 0)

Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4}; // Free Tip Face (X = L)

Curve Loop(5) = {1, 10, -5, -9};
Plane Surface(5) = {5}; // Bottom Face (Y = 0)

Curve Loop(6) = {3, 12, -7, -11};
Plane Surface(6) = {6}; // Top Face (Y = H)

// Define Volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Transfinite Mesh Setup for Hexahedral Elements
Transfinite Curve {1, 3, 5, 7} = 61; // 60 divisions along X
Transfinite Curve {2, 4, 6, 8} = 9;  // 8 divisions along Y
Transfinite Curve {9, 10, 11, 12} = 3; // 2 divisions along Z

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume(1);

// Physical Entities for Boundary Conditions
Physical Surface("Fixed_Support", 1) = {3}; // X = 0
Physical Surface("Tip_Face", 2) = {4};      // X = L
Physical Volume("Steel_Domain", 1) = {1};
