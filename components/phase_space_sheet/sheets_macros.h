#define SET_GAMMAI_DER(I) \
  d##I##gammai11_a(i, j, k) = -4.0*derivative(i, j, k, I, DIFFphi_a)*gammai11 \
+ std::exp(-4.0*DIFFphi_a(i, j, k))*(derivative(i, j, k, I, DIFFgamma22_a) + derivative(i, j, k, I, DIFFgamma33_a) - 2.0*derivative(i, j, k, I, DIFFgamma23_a) + derivative(i, j, k, I, DIFFgamma22_a)*DIFFgamma33_a(i, j, k) + DIFFgamma22_a(i, j, k)*derivative(i, j, k, I, DIFFgamma33_a)); \
    \
  d##I##gammai22_a(i, j, k) = -4.0*derivative(i, j, k, I, DIFFphi_a)*gammai22 \
  + std::exp(-4.0*DIFFphi_a(i, j, k))*(derivative(i, j, k, I, DIFFgamma11_a) + derivative(i, j, k, I, DIFFgamma33_a) - 2.0*derivative(i, j, k, I, DIFFgamma13_a) + derivative(i, j, k, I, DIFFgamma11_a)*DIFFgamma33_a(i, j, k) + DIFFgamma11_a(i, j, k)*derivative(i, j, k, I, DIFFgamma33_a));\
    \
   d##I##gammai33_a(i, j, k) = -4.0*derivative(i, j, k, I, DIFFphi_a)*gammai33 \
  + std::exp(-4.0*DIFFphi_a(i, j, k))*(derivative(i, j, k, I, DIFFgamma11_a) + derivative(i, j, k, I, DIFFgamma22_a) - 2.0*derivative(i, j, k, I, DIFFgamma12_a) + derivative(i, j, k, I, DIFFgamma11_a)*DIFFgamma22_a(i, j, k) + DIFFgamma11_a(i, j, k)*derivative(i, j, k, I, DIFFgamma22_a));\
    \
    d##I##gammai12_a(i, j, k) = -4.0*derivative(i, j, k, I, DIFFphi_a)*gammai12 \
  + std::exp(-4.0*DIFFphi_a(i, j, k))*(derivative(i, j, k, I, DIFFgamma13_a)*DIFFgamma23_a(i, j, k) + DIFFgamma13_a(i, j, k)*derivative(i, j, k, I, DIFFgamma23_a) - derivative(i, j, k, I, DIFFgamma12_a)*(1.0 + DIFFgamma33_a(i, j, k)) - DIFFgamma12_a(i, j, k)*derivative(i, j, k, I, DIFFgamma33_a)); \
    \
    d##I##gammai13_a(i, j, k) = -4.0*derivative(i, j, k, I, DIFFphi_a)*gammai13 \
  + std::exp(-4.0*DIFFphi_a(i, j, k))*(derivative(i, j, k, I, DIFFgamma12_a)*DIFFgamma23_a(i, j, k) + DIFFgamma12_a(i, j, k)*derivative(i, j, k, I, DIFFgamma23_a) - derivative(i, j, k, I, DIFFgamma13_a)*(1.0 + DIFFgamma22_a(i, j, k)) - DIFFgamma13_a(i, j, k)*derivative(i, j, k, I, DIFFgamma22_a)); \
\
    d##I##gammai23_a(i, j, k) = -4.0*derivative(i, j, k, I, DIFFphi_a)*gammai23 \
  + std::exp(-4.0*DIFFphi_a(i, j, k))*(derivative(i, j, k, I, DIFFgamma12_a)*DIFFgamma13_a(i, j, k) + DIFFgamma12_a(i, j, k)*derivative(i, j, k, I, DIFFgamma13_a) - derivative(i, j, k, I, DIFFgamma23_a)*(1.0 + DIFFgamma11_a(i, j, k)) - DIFFgamma23_a(i, j, k)*derivative(i, j, k, I, DIFFgamma11_a)) 

