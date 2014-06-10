#ifndef BSSN_MACROS
#define BSSN_MACROS

// applying functions to lots of vars

#define BSSN_APPLY_TO_FIELDS(function)  \
        function(gamma11);              \
        function(gamma12);              \
        function(gamma13);              \
        function(gamma22);              \
        function(gamma23);              \
        function(gamma33);              \
        function(gammai11);             \
        function(gammai12);             \
        function(gammai13);             \
        function(gammai22);             \
        function(gammai23);             \
        function(gammai33);             \
        function(phi);                  \
        function(A11);                  \
        function(A12);                  \
        function(A13);                  \
        function(A22);                  \
        function(A23);                  \
        function(A33);                  \
        function(K);                    \
        function(Gamma1);               \
        function(Gamma2);               \
        function(Gamma3);               \
        function(beta1);                \
        function(beta2);                \
        function(beta3);                \
        function(alpha);

#define BSSN_APPLY_TO_IJ_PERMS(function, more_args) \
        function(1, 1, more_args);                  \
        function(1, 2, more_args);                  \
        function(1, 3, more_args);                  \
        function(2, 2, more_args);                  \
        function(2, 3, more_args);                  \
        function(3, 3, more_args);

#define BSSN_APPLY_TO_JK_PERMS(more_args, function) \
        function(more_args, 1, 1);                  \
        function(more_args, 1, 2);                  \
        function(more_args, 1, 3);                  \
        function(more_args, 2, 2);                  \
        function(more_args, 2, 3);                  \
        function(more_args, 3, 3);


// BSSN & metric calculations

#define BSSN_CALCULATE_CHRISTOFFEL(I, J, K) G##I##J##K## = 0.5*( \
    gammai##I##1 * (d##J##g##K##1 + d##K##g##J##1 - d1g##J##K##) + \
    gammai##I##2 * (d##J##g##K##2 + d##K##g##J##2 - d2g##J##K##) + \
    gammai##I##3 * (d##J##g##K##3 + d##K##g##J##3 - d3g##J##K##) \
  )

#define BSSN_CALCULATE_DGAMMAI(I, J, K) d##I##gi##J##K## = der(gammai##J##K##, I, idx);

#define BSSN_CALCULATE_DGAMMA(I, J, K) d##I##g##J##K## = der(gamma##J##K##, I, idx);

// standard ordering of indexes for tensor components

// actual fields:
#define gamma21 gamma12
#define gamma31 gamma13
#define gamma32 gamma23
#define gammai21 gammai12
#define gammai31 gammai13
#define gammai32 gammai23
#define A21 A12
#define A31 A13
#define A32 A23

// local variables:
// ricci tensor
#define ricciTF21 ricciTF12
#define ricciTF31 ricciTF13
#define ricciTF32 ricciTF23

// covariant double-derivatives of alpha
#define D2D1a D1D2a
#define D3D1a D1D3a
#define D3D2a D2D3a

// covariant double-derivatives of phi
#define D2D1phi D1D2phi
#define D3D1phi D1D3phi
#define D3D2phi D2D3phi

// christoffel symbols
#define G121 G112
#define G131 G113
#define G132 G123
#define G221 G212
#define G231 G213
#define G232 G223
#define G321 G312
#define G331 G313
#define G332 G323

// Metric derivatives
#define d1g21 d1g12
#define d1g31 d1g13
#define d1g32 d1g23
#define d2g21 d2g12
#define d2g31 d2g13
#define d2g32 d2g23
#define d3g21 d3g12
#define d3g31 d3g13
#define d3g32 d3g23

// inverse metric derivatives
#define d1gi21 d1gi12
#define d1gi31 d1gi13
#define d1gi32 d1gi23
#define d2gi21 d2gi12
#define d2gi31 d2gi13
#define d2gi32 d2gi23
#define d3gi21 d3gi12
#define d3gi31 d3gi13
#define d3gi32 d3gi23

// second derivatives of the metric
// bad metric indices
#define d1d1g21 d1d1g12
#define d1d1g31 d1d1g13
#define d1d1g32 d1d1g23
#define d1d2g21 d1d2g12
#define d1d2g31 d1d2g13
#define d1d2g32 d1d2g23
#define d1d3g21 d1d3g12
#define d1d3g31 d1d3g13
#define d1d3g32 d1d3g23
#define d2d2g21 d2d2g12
#define d2d2g31 d2d2g13
#define d2d2g32 d2d2g23
#define d2d3g21 d2d3g12
#define d2d3g31 d2d3g13
#define d2d3g32 d2d3g23
#define d3d3g21 d3d3g12
#define d3d3g31 d3d3g13
#define d3d3g32 d3d3g23
// bad derivative indices
#define d2d1g11 d1d2g11
#define d2d1g12 d1d2g12
#define d2d1g13 d1d2g13
#define d2d1g22 d1d2g22
#define d2d1g23 d1d2g23
#define d2d1g33 d1d2g33
#define d3d1g11 d1d3g11
#define d3d1g12 d1d3g12
#define d3d1g13 d1d3g13
#define d3d1g22 d1d3g22
#define d3d1g23 d1d3g23
#define d3d1g33 d1d3g33
#define d3d2g11 d2d3g11
#define d3d2g12 d2d3g12
#define d3d2g13 d2d3g13
#define d3d2g22 d2d3g22
#define d3d2g23 d2d3g23
#define d3d2g33 d2d3g33
// bad in both indices
#define d2d1g21 d1d2g12
#define d2d1g31 d1d2g13
#define d2d1g32 d1d2g23
#define d3d1g21 d1d3g12
#define d3d1g31 d1d3g13
#define d3d1g32 d1d3g23
#define d3d2g21 d2d3g12
#define d3d2g31 d2d3g13
#define d3d2g32 d2d3g23

#endif
