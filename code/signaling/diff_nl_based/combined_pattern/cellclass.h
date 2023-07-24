// Cell class define
#ifndef _CELL_H_
#define _CELL_H_
//#include <iostream>
//#include <fstream>
#include <string>
#include <iostream>

//using std::string;

#include "decvar.h"
#include <vector>

using namespace std;

struct GCpara
{
  int N;                   // number of cells
  int Nc;                 // overall maximum number of c-sites formed between two cells
  double r;                // radius of the cell
  double Lc;               // maximum length of the spring for c-sites 
  vector<double> argC;     // direction of cAMP at each cell
  vector<double> rangeC;           // reach-out range
  double tc_PP;            // attach time for c-sites
  double tc_TT;
  double tc_TP;
};

struct GIpara
{
  int N;                   // number of cells
  int Ni;                  // maximum number of c-sites formed between two cells
  double r;                // radius of the cell
  double Li;               // maximum length of the spring for c-sites 
  vector<double> argC;     // direction of cAMP
  double ti_PST;           // attach time for i-sites of PST celss
  double ti_PSP;           // attach time for i-sites of PSP celss
};

struct Mpara
{
  int N;                                 // number of cells
  int Ni;                                // maximum number of c-sites formed between two cells
  double r;                              // radius of the cell
  double lc;                             // natural length of the spring for c-sites 
  double Lc;                             // maximum length of the spring for c-sites 
  double nl;                             // differential natural length
  double li;                             // natural length of the spring for i-sites
  double Li;                             // maximum length of the spring for i-sites 
  vector<double> AlphaC;                 // spring constant of c-sites TT,PP,TP
  vector<double> AlphaI;                 // spring constant of i-sites TS,PS
  double alphaB;                         // spring constant for cell cell body force
  double alphaS;                         // spring constant for cell substrate body force
  vector<double> AlphaM;                 // linear body force on cell membrane
  double mu;                             // viscosity of cells
  double nu;                             // viscosity of c-sites: nu = mu/ratio
  vector<vector<int> > nc;               // pairwise number of c-sites 
  vector<int> ni;                        // number of i-sites of each cell
  vector<vector<vector<double> > > Si;   // all the i-sites locations
  vector<int> label;                     // label (PST/PSP) 0/1 of each cell
};

struct PDEpara
{
  int N;                            // number of cells
  double Lx;                         // length of periodic domain
  double Ly;                         // height of periodic domain
  int Nx; int Ny;                   // number of subintervals
  double mindt;                     // smallest timestep in CK algorithm
  double sigma;                     // coeff of the diffusion term
  double tol;                       // tolerance in CK algorithm
  double baseT; double baseP;       // baseal level of outputting chemicals
  double tshT; double tshP;         // threshold to output chemicals
  double Tt; double Tp;             // period of chemical output
  double slopT; double slopP;       // slope of chemical output
};

struct NLpara
{
  int N;                                 // number of cells
  int Nx; int Ny;                        // number of subintervals
  int K;                                 // weight factor
  double Lx;                              // length of periodic domain
  double Ly;                         // height of periodic domain
  double Vc;                             // volume of a cell
  double V0;                             // volume of the extracellular medium
  double src; vector<int> src_loc;    // exogenous source and source location
  double gm6_t; double gm6_p;
  double gm7_t; double gm7_p;            // local degradation
  double gm8; double gm9;                // global degradation
};

class Cadherin
{
 public:
	void gen_c(vector<vector<double> >& X, GCpara gc_para, vector<vector<vector<vector<double> > > >& Sc, vector<vector<int> >& nc, vector<vector<vector<double> > >& Tc, vector<int>& LB);
  void gen_c_init(vector<vector<double> >& X, GCpara gc_para, vector<vector<vector<vector<double> > > >& Sc, vector<vector<int> >& nc, vector<vector<vector<double> > >& Tc, vector<int>& LB);
  void gen_c_single(vector<vector<double> >& X, int i, int j, GCpara gc_para, double* Sc, double& Tc, vector<int>& LB);
  void removeC(vector<vector<vector<vector<double> > > >& Sc, vector<vector<int> >& nc, vector<vector<vector<double> > >& Tc, int* idx, int N, int Nc);
  void addC(vector<vector<vector<vector<double> > > >& Sc, vector<vector<int> >& nc, vector<vector<vector<double> > >& Tc, int* cidx, double* ccor, double ctime, int N, int Nc);
};

class Integrin
{
  public:
   void gen_i(vector<vector<double> >& X, GIpara gi_para, vector<vector<vector<double> > >& Si, vector<int>& ni, vector<vector<double> >& Ti, vector<int>& LB);
   void gen_i_init(vector<vector<double> >& X, GIpara gi_para, vector<vector<vector<double> > >& Si, vector<int>& ni, vector<vector<double> >& Ti, vector<int>& LB);
   void gen_i_single(vector<vector<double> >& X, int tba_i, GIpara gi_para, double* Si, double& Ti, vector<int>& LB);
   void removeI(vector<vector<vector<double> > >& Si, vector<int>& ni, vector<vector<double> >& Ti, int* idx, int N, int Ni);
   void addI(vector<vector<vector<double> > >& Si, vector<int>& ni, vector<vector<double> >& Ti, int idx, double* icor, double itime, int N, int Ni);
};


#endif
