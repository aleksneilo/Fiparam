#include <complex>

 using namespace std;
 
 extern int N_SL, N_SR, N_Mid, N, N_S, NUM_Tech, intG, idGm, idGm1;
 extern double epsG, epsDel, alpha; 
extern int  w_obrez, iter;
extern int  MODE;  // regime of boundary conditions: 0 - Free Boundary, 1 - Fixed phase
extern complex<double> *difG, *difGm;
 extern double T, Del0, Del0a, Del0b,Xi1, Xi2;
 extern double L_SL, L_SR, L_S, L_F1, Ksi_S, ro_S, Tc;
 extern double L_F, Ksi_F, ro_F, H, alphaT, aGIP, aGIPlast;
 extern double Gb_I, Gb_FS, Gb_SF, R0A;
 extern double *Hi, *Ksii, *Roi, *Li, *Stype, *Tci, *Rbi, *a1coef, k1, k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14;
 
 
 void SelfCons(complex<double> *G, complex<double> *Del, int Initial);
 double SelfConsZero();
 void getABC(complex<double> *a, complex<double> *b, complex<double> *c, complex<double> *d, int iN, complex<double> *G, complex<double> *Del, complex<double> w);
 void Prog( complex<double> *Fi, complex<double> *G, complex<double> *Del, complex<double> w);
 void Gcalc( complex<double> *G, double *dGmax, complex<double> *Fi1, complex<double> *Fi2, complex<double> *Del, complex<double> w);
 void GcalcDOS(complex<double> *G, double *dGmax, complex<double> *Fi1, complex<double> *Fi2, complex<double> *Del, complex<double> w);
 double E_IS_calc(double *Is, complex<double> *G,complex<double> *Del);
 void DOS(double E, complex<double> *G, complex<double> *G1, complex<double> *Fi, complex<double> *Fi1, complex<double> *Del);
 complex<double> Gcalc0(int i, complex<double> w, complex<double> *Del);
 void Discrepancy( complex<double> *Fi, complex<double> *G, complex<double> *Del, complex<double> w, double *rmax);
 complex<double> SQRT(int i, complex<double> EE, complex<double> ehtr, complex<double> *Del);
 double sign(double x);
 
 int Layer(int iN);
 complex<double> get_wm(int iN, complex<double> w);
 double get_h(int iN);
 double get_HE(int iN);
 double get_h_out(int iN);
 double get_ksi(int iN);
 double get_type(int iN);
 double get_tc(int iN);
 

