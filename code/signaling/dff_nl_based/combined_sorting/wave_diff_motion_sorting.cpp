#include <fstream>
#include <cmath>
#include <random>
#include <iterator>
#include <string>
#include "cellclass.h"
#include <cvode/cvode.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <nvector/nvector_serial.h>
#include <vector>
// #include <lapack.h>
#include <Accelerate/Accelerate.h>
// #include <clapack.h>
#include <sundials/sundials_types.h>
using namespace std;

/************************************************************************************
                           function declarations
************************************************************************************/
int check_flag(void *flagvalue, char *funcname, int opt);
double norm(double *x, int dim);
double sigmoid(double x, double m, double c, double s);
void redim21(double (*V)[2], int len_V, double *v);
void redim12(double *v, int len_v, double (*V)[2]);
int pre_count(int i, int j, vector<vector<int>> &nc, int N);
void gen_init(int N, vector<vector<double>> &X, vector<vector<int>> &nc, vector<vector<vector<vector<double>>>> &Sc, double *cc_array);
void Id_c(int nij, double *xi, double *xj, vector<vector<double>> &sij, double r, double L, vector<double> &Phi);
void f_red(int nij, double alpha, double r, double lc, double Lc, double nl, double *x, vector<vector<double>> &s, vector<vector<double>> &F);
void f_blue(int nij, double alpha, double r, double lc, double Lc, double *x, vector<vector<double>> &s, vector<vector<double>> &F);
void f1_red(double alpha, double r, double lc, double Lc, double nl, double *x, double *y, double *F);
void f1_blue(double alpha, double r, double lc, double Lc, double *x, double *y, double *F);
double Id_b(double *xi, double *xj, double r, double lc);
void g(double alpha, double alphaM, double r, double lc, double *x, double *y, double *B);
int f_rhs(double t, N_Vector y, N_Vector ydot, void *user_data);
void ADI_nonlin(vector<double> &u_init, vector<double> &sr_init, double tstart, double tend, vector<vector<double>> &X, vector<int> label, vector<double> &tchem, PDEpara CKpara, NLpara Fpara, vector<double> &u, vector<double> &rtchem, vector<double> &sr_out, vector<int> &pulse_marker);
void nonlin_F(vector<double> &upre, vector<double> &sr, vector<vector<double>> &X, vector<int> &label, NLpara para, vector<double> &x, vector<double> &y, vector<double> &F);
// random sampling functions
template <typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator &g)
{
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}
template <typename Iter>
Iter select_randomly(Iter start, Iter end)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}
/************************************************************************************
                                global variables
************************************************************************************/
double tol = 1e-4;       // tolerance for steady states
double T_MAX = 180;      // maximum time for motion equations
double r = 1;            // radius of the nucleus (micrometer)
double lc = 2;           // natural length of the spring for c-sites
double li = 1;           // natural length of the spring for i-sites
double Lc = 8;           // max length of the spring for c-sites (micrometer)
double nl = 2;           // diff natural length
double Li = 7;           // max length of the spring for i-sites (micrometer)
double argC = 0;         // angle of the gradient of cAMP concentration [-pi,pi]
double rangeC = pi / 16; // reach out range
int N = 48;              // total number of cells
int rowN = 12;           // number of cells in each row
int columnN = 4;         // number of cells in each column
double gap = 2 * r;
double margin_x = 10 * r;
double margin_y = 0 * r;
int Nc = 4;                           // max number of c-sites formed between two cells
int mNc = 2;                          // min number of c-sites formed between two cells
int Ni = 4;                           // max number of i-sites one cell can form
int mNi = 1;                          // min number of i-sites of a cell close to the substrate
vector<double> AlphaC{350, 350};      // spring constant of c-sites (g/min^2), T, P
vector<double> AlphaI{0, 350};        // spring constant of i-sites (g/min^2)
double alphaB = 2000;                 // spring constant for cell cell body force (g/min^2)
vector<double> AlphaM = {50, 50, 50}; // linear body force on cell membrane, TT, PP, TP
double alphaS = 2500;                 // spring constant for cell substrate body force (g/min^2)
double vis_ratio = 10;
double mu = 0.101 / 60;     // viscosity of cells (g/min)
double nu = mu / vis_ratio; // viscosity of c-sites (g/min)
double tc_TT = 3;           // attach time for each c-site (min)
double tc_PP = 3;
double tc_TP = 3;
double ti_PSP = 0.1; // attach time for each i-site (min)
double ti_PST = 0.1;

// constants for ODE solver
const int kry_maxl = 4; // dimension of krylov space
double reltol = 1e-6;
double abstol = 1e-6;
int maxstep = 10000; // maximum steps to be taken for the ode solver

// discretize spatial domain for pde solver
double Lx = margin_x + 2 * r * rowN + margin_x; // length of the periodic domain
double Ly = margin_y + 2 * r * columnN + 4;     // height of the periodic domain
int Nx = 2 * Lx;                                // number of x-subintervals
int Ny = 2 * Ly;                                // number of y-subintervals
double chaL = 5.5 * r;                          // characteristic dimension of the system
int mx = Nx + 1;
int my = Ny + 1;
double dx = Lx / Nx;
double dy = Ly / Ny;                                      // length of subintervals
static const long long int MAX_ENTRY = 89 * 25 * 89 * 25; // # of matrix entries in the linear system
static const int MAX_DIM = 89 * 25;                       // dimension of rhs in the linear equation
static const int NX = 88;
static const int NY = 24;
double A_longarray[MAX_ENTRY] = {0};
double rhs_longarray[MAX_DIM] = {0};
float Dx[NX - 1] = {0};
float Dy[NY - 1] = {0};
float Bx[NX - 2] = {0};
float By[NY - 2] = {0};

// constants for customized PDE solver
double mindt = 1e-5;                  // smallest timestep in CK algorithm
double CKtol = 1e-3;                  // tolerance in CK algorithm
double sigma_c = 300 / (chaL * chaL); // diffusion coeff of camp
double sigma_d = 300 / (chaL * chaL); // diffusion coeff of dif
double iPDE = 1.73;                   // [iPDE]_T (micromolars)
double T_ct = 2;
double T_cp = 2; // period of camp output
double T_dt = 2;
double T_dp = 2;               // period of dif output
double base_ct = 0.099 / iPDE; // modified                    // basal level of PST outputting camp
double base_cp = base_ct / 10; // basal level of PSP outputtting camp
double base_dt = 0;            // basal level of PST outputting dif
double base_dp = base_ct;      // basal level of PSP outputtting dif
double slop_ct = 0;            // 47.2 / iPDE;                           // slope of PST outputting camp (/min)
double slop_cp = slop_ct / 10; // slope of PSP outputting camp (/min)
double slop_dt = 0;            // slope of PST outputting dif (/min)
double slop_dp = 0;            // slope of PSP outputting dif (/min)
double tsh_ct = 5 * 1e-2;      // derivative threshold for PST outputting camp
double tsh_cp = 5 * 1e-2;      // derivative threshold for PSP outputting camp
double tsh_dt = 5 * 1e-3;      // derivative threshold for PST outputting dif
double tsh_dp = 5 * 1e-3;      // derivative threshold for PSP outputting dif

// constants for nonlinear term in pde
double src_c = base_ct * 10; // constant exogenous camp source
double src_d = 0;            // constant exogenous dif source
// int src_loc_c[2] = [Nx - Nx/5, Ny - Ny/5];            // cAMP source location
// vector<int> src_loc_d{1,1};                             // dif source location
int K_w = 1;                                  // weight factor
double V0 = (Lx * chaL) * (Ly * chaL) * chaL; // modified                             // volume of the extracellular medium
double Vc = 696.9;                            // volume of a cell (micrometers^3)
double gm6_ct = 11.6;
double gm6_cp = 11.6;
double gm6_dt = 1.16;
double gm6_dp = 1.16;
double gm7_ct = 0;
double gm7_cp = 0; // no local degradation
double gm7_dt = gm7_cp;
double gm7_dp = 0; // no local degradation
double gm8_c = 0.1;
double gm8_d = gm8_c;
double gm9_c = 0.1;
double gm9_d = gm9_c; // small global degradation
double gm9_hat_c = gm9_c * N * K_w * Vc / V0;
double gm9_hat_d = gm9_hat_c; // modified

// parameters for cell differentiation
double tsh_PT = 0.25; // feedback threshold for PSP turning into PST
double tsh_TP = 0.25; // feedback threshold for PST turning into PSP
double a0 = 0.0005;   // half-max dif for PSP -> PST (micromolars)
double b0 = 0.0002;   // half-max dif for PST -> PSP (micromolars)
double a1 = 0.0004;   // half-max camp for PST -> PSP (micromolars)
double b1 = 0.0002;   // half-max camp for PSP -> PST (micromolars)
double exp_c = 1;
double exp_d = 1;     // exponents in transition functions
double int_fdb_t = 1; // start collecting feedback info 1 minute prior to differentiation time
vector<double> Feedback_tp(N, 0);
vector<double> Feedback_pt(N, 0);

int main()
{
    std::random_device rd;
    std::mt19937 generator(rd());

    SUNContext sunctx;
    int sunflag = 0;
    sunflag = SUNContext_Create(NULL, &sunctx);

    // call structures & classes
    GCpara gc_para;
    GIpara gi_para;
    Mpara motion_para;
    Mpara *mpara = &motion_para;
    PDEpara PCpara;
    PDEpara PDpara;
    NLpara FCpara;
    NLpara FDpara;
    Cadherin cad;
    Integrin intg;
    SUNLinearSolver LS;

    // parameters for generating c-sites
    gc_para.N = N;
    gc_para.Nc = Nc;
    gc_para.r = r;
    gc_para.Lc = Lc;
    gc_para.tc_PP = tc_PP;
    gc_para.tc_TT = tc_TT;
    gc_para.tc_TP = tc_TP;

    // parameters for motion equations
    motion_para.N = N;
    motion_para.Ni = Ni;
    motion_para.r = r;
    motion_para.lc = lc;
    motion_para.Lc = Lc;
    motion_para.nl = nl;
    motion_para.AlphaC = AlphaC;
    motion_para.alphaB = alphaB;
    motion_para.alphaS = alphaS;
    motion_para.AlphaM = AlphaM;
    motion_para.mu = mu;
    motion_para.nu = nu;

    // parameters for cAMP pde solver
    PCpara.N = N;
    PCpara.Lx = Lx;
    PCpara.Ly = Ly;
    PCpara.Nx = Nx;
    PCpara.Ny = Ny;
    PCpara.mindt = mindt;
    PCpara.sigma = sigma_c;
    PCpara.tol = CKtol;
    PCpara.baseT = base_ct;
    PCpara.baseP = base_cp;
    PCpara.tshT = tsh_ct;
    PCpara.tshP = tsh_cp;
    PCpara.Tt = T_ct;
    PCpara.Tp = T_cp;
    PCpara.slopT = slop_ct;
    PCpara.slopP = slop_cp;

    // parameters for DIF pde solver
    PDpara.N = N;
    PDpara.Lx = Lx;
    PDpara.Ly = Ly;
    PDpara.Nx = Nx;
    PDpara.Ny = Ny;
    PDpara.mindt = mindt;
    PDpara.sigma = sigma_d;
    PDpara.tol = CKtol;
    PDpara.baseT = base_dt;
    PDpara.baseP = base_dp;
    PDpara.tshT = tsh_dt;
    PDpara.tshP = tsh_dp;
    PDpara.Tt = T_dt;
    PDpara.Tp = T_dp;
    PDpara.slopT = slop_dt;
    PDpara.slopP = slop_dp;

    // parameters for nonlinear term in cAMP pde
    FCpara.N = N;
    FCpara.Nx = Nx;
    FCpara.Ny = Ny;
    FCpara.K = K_w;
    FCpara.Lx = Lx;
    FCpara.Ly = Ly;
    FCpara.Vc = Vc;
    FCpara.V0 = V0;
    // FCpara.src = src_c; // FCpara.src_loc = src_loc_c;
    FCpara.gm6_t = gm6_ct;
    FCpara.gm6_p = gm6_cp;
    FCpara.gm7_t = gm7_ct;
    FCpara.gm7_p = gm7_cp;
    FCpara.gm8 = gm8_c;
    FCpara.gm9 = gm9_hat_c;

    // parameters for nonlinear term in DIF pde
    FDpara.N = N;
    FDpara.Nx = Nx;
    FDpara.Ny = Ny;
    FDpara.K = K_w;
    FDpara.Lx = Lx;
    FDpara.Ly = Ly;
    FDpara.Vc = Vc;
    FDpara.V0 = V0;
    // FDpara.src = src_d; FDpara.src_loc = src_loc_d;
    FDpara.gm6_t = gm6_dt;
    FDpara.gm6_p = gm6_dp;
    FDpara.gm7_t = gm7_dt;
    FDpara.gm7_p = gm7_dp;
    FDpara.gm8 = gm8_d;
    FDpara.gm9 = gm9_hat_d;

    // grid points including boundaries
    vector<double> x(mx, 0);
    vector<double> y(my, 0);
    for (int j = 0; j < mx; j++)
    {
        x[j] = j * dx;
    }
    for (int k = 0; k < my; k++)
    {
        y[k] = k * dy;
    }

    /*****************************************************************
                                Time Tables
    ******************************************************************/
    // time table for the system to generate new cadherins
    vector<double> Tgc = {};
    vector<double> Tgi = {};
    double cumultimeC = 0;
    double cumultimeI = 0;
    int Tgc_len = 0;
    int Tgi_len = 0;

    while (cumultimeC < T_MAX)
    {
        std::normal_distribution<> NdistrC(0.01, 0.005);
        double gc_interval = NdistrC(generator);
        if (gc_interval > 0)
        {
            cumultimeC += gc_interval;
            Tgc.push_back(cumultimeC);
            Tgc_len += 1;
        }
    }

    // check cell type every eight minutes (time table for cell differentiation)
    vector<int> Tlb(round(T_MAX) / 8 + 1, 0);
    for (int i = 0; i <= round(T_MAX) / 8; i++)
    {
        Tlb[i] = 20 + i * 8; // freeze for 20 minutes
    }
    int lb_check = 0;

    /*****************************************************************
                                Initialization
    ******************************************************************/
    // create output files ------------------------------------------
    string filename1 = "motion.dat";
    ofstream *foutm = new ofstream(filename1);
    string filename2 = "dif.dat";
    ofstream *foutd = new ofstream(filename2);
    string filename3 = "camp.dat";
    ofstream *foutc = new ofstream(filename3);

    // initialize cell centers --------------------------------------
    vector<vector<double>> X(N, vector<double>(2, 0));
    // initial rectangular shape
    for (int i = 0; i < N; i++)
    {
        int posx = fmod(i, rowN);
        int posy = i / rowN;
        X[i][0] = r + margin_x + posx * gap;
        X[i][1] = r + margin_y + posy * gap;
    }

    // initialize cell types -----------------------------------------
    vector<int> LB(N, 1); // 0: PST; 1: PSP
    for (int i = 0; i < 4; i++)
    {
        int lbpos1 = 12 * i;
        int lbpos2 = lbpos1 + 1;
        int lbpos3 = lbpos2 + 1;
        LB[lbpos1] = 0;
        LB[lbpos2] = 0;
        LB[lbpos3] = 0;
    }
    lb_check += 1;

    // initialize signalling settings --------------------------------
    vector<double> uc(mx * my, 0);
    vector<double> ud(mx * my, 0);
    vector<double> sr_init_c(N, 0);
    vector<double> sr_init_d(N, 0);
    vector<double> sr_out_c(N, 0);
    vector<double> sr_out_d(N, 0);

    vector<int> pulse_marker_c(N, 0); // initialize pulsing state
    int src_c_idx = 35;               // int src_d_idx = 24; // fix source cell

    for (int i = 0; i < N; i++)
    {
        if (LB[i] == 0)
        {
            sr_init_c[i] = base_ct; // base_ct
        }
        else
        {
            sr_init_c[i] = base_cp; // base_cp
        }
    }
    sr_out_c = sr_init_c;
    vector<double> rtcamp(N, 0);
    vector<double> rtdif(N, 0);

    // write initial camp & dif concentration
    double center = margin_x + 9 * gap; // x = 30
    double scale = r * 8 / 8;
    vector<vector<double>> ruc(mx, vector<double>(my, 0));
    for (int j = 0; j < mx; j++)
    {
        for (int k = 0; k < my; k++)
        {
            int pos = k * mx + j;
            if (x[j] >= X[0][0] && x[j] <= X[47][0] && y[k] >= X[0][1] && y[k] <= X[47][1])
            {
                ud[pos] = sigmoid(x[j], base_dp, center, scale);
                uc[pos] = sigmoid(-x[j], base_ct, -center, scale);
            }
            ruc[j][k] = uc[pos];
            *foutd << ud[pos] << " ";
            *foutc << uc[pos] << " ";
        }
        *foutd << "\n";
        *foutc << "\n";
    }

    // initial argC direction of camp -----------------------------------
    vector<double> ArgC(N, argC);
    vector<double> ux(N, 0);
    vector<double> uy(N, 0); // partial derivative
    vector<vector<double>> Xmod(N, vector<double>(2, 0));
    for (int i = 0; i < N; i++)
    {
        Xmod[i][0] = std::fmod(X[i][0], Lx);
        Xmod[i][1] = std::fmod(X[i][1], Ly);
    }
    for (int i = 0; i < N; i++)
    {
        double cellP[2] = {0, 0};
        cellP[0] = Xmod[i][0];
        cellP[1] = Xmod[i][1]; // the i-th cell's location
        double dist = Lx;
        int k_clo = 0;
        int j_clo = 0;
        for (int k = 0; k < my; k++)
        {
            for (int j = 0; j < mx; j++)
            {
                double gridP[2] = {x[j], y[k]};
                double tmp_dist = 0;
                for (int l = 0; l < 2; l++)
                {
                    tmp_dist += (cellP[l] - gridP[l]) * (cellP[l] - gridP[l]);
                }
                tmp_dist = sqrt(tmp_dist);
                if (tmp_dist < dist)
                {
                    dist = tmp_dist;
                    k_clo = k;
                    j_clo = j;
                }
            }
        } //  find the grid point closest to cell i
        if (j_clo == 0)
        {
            ux[i] = (ruc[j_clo + 1][k_clo] - ruc[j_clo][k_clo]) / dx;
        }
        else if (j_clo == Nx)
        {
            ux[i] = (ruc[j_clo][k_clo] - ruc[j_clo - 1][k_clo]) / dx;
        }
        if (k_clo == 0)
        {
            uy[i] = (ruc[j_clo][k_clo + 1] - ruc[j_clo][k_clo]) / dy;
        }
        else if (k_clo == Ny)
        {
            uy[i] = (ruc[j_clo][k_clo] - ruc[j_clo][k_clo - 1]) / dy;
        }
        if (j_clo != 0 && k_clo != 0 && j_clo != Nx && k_clo != Ny)
        {
            ux[i] = (ruc[j_clo + 1][k_clo] - ruc[j_clo - 1][k_clo]) / (2 * dx);
            uy[i] = (ruc[j_clo][k_clo + 1] - ruc[j_clo][k_clo - 1]) / (2 * dy);
        }
        ArgC[i] = atan2(uy[i], ux[i]); // lie in [-pi,pi]
    }
    gc_para.argC = ArgC;

    // initialize all the cadherins -----------------------------------
    vector<vector<int>> nc(N, vector<int>(N, 0));                                                                                      // int n_c[N * N];
    vector<vector<vector<vector<double>>>> Sc(N, vector<vector<vector<double>>>(N, vector<vector<double>>(Nc, vector<double>(2, 0)))); // double S_c[Tnc][2];
    vector<vector<vector<double>>> Tc(N, vector<vector<double>>(N, vector<double>(Nc, 0)));                                            // double T_c[Tnc];

    vector<double> RangeC(N, rangeC);
    gc_para.rangeC = RangeC;

    cad.gen_c_init(X, gc_para, Sc, nc, Tc, LB);

    // write initialization ----------------------------------
    // write initial nc
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            *foutm << nc[i][j] << " ";
        }
        *foutm << "\n";
    }

    // wirte initial cell labels
    for (int i = 0; i < N; i++)
    {
        *foutm << LB[i] << " ";
    }
    *foutm << "\n";

    // write initial cells and c-sites locations
    int Tnc_init = pre_count(N - 1, N, nc, N);
    N_Vector init_cond_init;
    double init_cond_array_init[10000];
    gen_init(N, X, nc, Sc, init_cond_array_init);
    init_cond_init = N_VMake_Serial((N + Tnc_init) * 2, init_cond_array_init, sunctx);
    for (int i = 0; i < (N + Tnc_init) * 2; i++)
    {
        if (i % 2 == 0)
        {
            *foutm << NV_Ith_S(init_cond_init, i) << " " << NV_Ith_S(init_cond_init, i + 1) << "\n";
        }
    }

    /************************************************************************************
                                    Stochastic System
    ************************************************************************************/
    int flag;
    double time_count = 0;

    while (time_count < T_MAX)
    {
        // create CVode object
        void *cvode_mem;
        cvode_mem = NULL;

        SUNLinearSolver LS;
        cvode_mem = CVodeCreate(CV_BDF, sunctx);

        // initial conditions
        int Tnc = pre_count(N - 1, N, nc, N);
        N_Vector init_cond;
        double init_cond_array[10000];
        gen_init(N, X, nc, Sc, init_cond_array);
        init_cond = N_VMake_Serial((N + Tnc) * 2, init_cond_array, sunctx);

        // set pacemaker, pulse every three minutes
        if (fmod(time_count, 4) >= 0 && fmod(time_count, 4) <= 1)
        {
            src_c = base_ct * 10;
            FCpara.src = src_c;
        }
        else
        {
            src_c = 0;
            FCpara.src = src_c;
        }

        // the cadherin that wears out the fastest
        double t1 = 100;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < Nc; k++)
                {
                    if (Tc[i][j][k] > 0 && Tc[i][j][k] < t1)
                    {
                        t1 = Tc[i][j][k];
                    }
                }
            }
        }
        if (t1 == 100)
        { // if no entry in Tc is positive
            t1 = 0;
            cout << "no entry in Tc is positive"
                 << "\n";
        }
        cout << "t1 = " << t1 << "\n";

        // the most recent attempt to form a new c-site
        vector<double> Tgc_left(Tgc_len, 0);
        for (int i = 0; i < Tgc_len; i++)
        {
            Tgc_left[i] = Tgc[i] - time_count;
        }
        double t2 = 100;
        for (int i = 0; i < Tgc_len; i++)
        {
            if (Tgc_left[i] > 0 && Tgc_left[i] < t2)
            {
                t2 = Tgc_left[i];
            }
        }
        if (t2 == 100)
        { // if no entry in Tgc_left is positive
            t2 = 0;
            cout << "no entry in Tgc_left is positive"
                 << "\n";
        }
        cout << "t2 = " << t2 << "\n";

        // determine the next step: remove/add a c-site/i-sites
        double t0 = 100;
        if (t1 > 0 && t1 < t0)
        {
            t0 = t1;
        }
        if (t2 > 0 && t2 < t0)
        {
            t0 = t2;
        }

        if (t0 == 100)
        {
            // no more attempt to form a new c/i-site
            // or no attached c/i-site left
            cout << "t0 = 0"
                 << "\n";
            break;
        }
        string action;
        if (t0 == t1)
        {
            action = "rc";
        }
        else if (t0 == t2)
        {
            action = "ac";
        }

        // cell motion during one time period for a fixed group of Sc
        motion_para.nc = nc;
        motion_para.label = LB;

        double timestep = t0 / 100;
        cout << "timestep: " << timestep << "\n";

        // write nc for this generation
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                *foutm << nc[i][j] << " ";
            }
            *foutm << "\n";
        }

        // write Tc for this generation
        for (int i = 0; i < N; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                for (int k = 0; k < Nc; k++)
                {
                    *foutm << Tc[i][j][k] << " ";
                }
                *foutm << "\n";
            }
        }

        // write how long this fixed group of c-sites & i-sites last
        *foutm << t0 << "\n";

        // wirte cell labels
        for (int i = 0; i < N; i++)
        {
            *foutm << LB[i] << " ";
        }
        *foutm << "\n";

        // write pulse markers
        cout << "pulse marker:"
             << " ";
        for (int i = 0; i < N; i++)
        {
            *foutm << pulse_marker_c[i] << " ";
            cout << pulse_marker_c[i] << " ";
        }
        *foutm << "\n";
        cout << "\n";

        if (time_count >= 20)
        {
            // save cell locations from each loop
            // vector<vector<vector<double> > > CellMove(N, vector<vector<double> >(2, vector<double>(100, 0)));
            vector<vector<double>> Xmod(N, vector<double>(2, 0));
            vector<vector<double>> rud(mx, vector<double>(my, 0));
            vector<vector<double>> ruc(mx, vector<double>(my, 0));
            vector<int> triggerState(N, 0);

            for (int ti = 0; ti < 100; ti++)
            {

                /*****************************************************************
                                        Motion Equations
                                        adcance by t0
                ******************************************************************/
                double tstart = time_count + ti * timestep;

                // reinitialize CVODE solver
                if (ti == 0)
                {
                    flag = CVodeInit(cvode_mem, f_rhs, tstart, init_cond); // the very first initialization
                }
                else
                {
                    flag = CVodeReInit(cvode_mem, tstart, init_cond);
                }
                CVodeSetUserData(cvode_mem, mpara);

                LS = SUNLinSol_SPGMR(init_cond, PREC_NONE, kry_maxl, sunctx);
                flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);

                flag = CVodeSetMaxNumSteps(cvode_mem, maxstep);
                flag = CVodeSStolerances(cvode_mem, reltol, abstol);

                // advance solution in time
                double tend = tstart + timestep;
                double tret;
                flag = CVode(cvode_mem, tend, init_cond, &tret, CV_ONE_STEP);

                // extract current cells locations
                for (int i = 0; i < N; i++)
                {
                    X[i][0] = NV_Ith_S(init_cond, 2 * i);
                    X[i][1] = NV_Ith_S(init_cond, 2 * i + 1);
                }
                for (int i = 0; i < N; i++)
                {
                    Xmod[i][0] = std::fmod(X[i][0], Lx);
                    Xmod[i][1] = std::fmod(X[i][1], Ly);
                }

                /*****************************************************************
                                    Signalling Equations
                        happen simultaneously with motion equations
                ******************************************************************/

                // dif equations --------------------------------------------
                vector<double> tdif(N, 0);
                tdif = rtdif;
                vector<double> ud_init(mx * my, 0);
                ud_init = ud;
                vector<int> pulse_marker_d(N, 0);
                // save preceeding result in ud_init and tdif
                // save output in ud and rtdif
                vector<double> sr_init_d(N, 0);
                sr_init_d = sr_out_d;
                ADI_nonlin(ud_init, sr_init_d, tstart, tend, X, LB, tdif, PDpara, FDpara, ud, rtdif, sr_out_d, pulse_marker_d);

                // camp equations -------------------------------------------
                vector<double> tcamp(N, 0);
                tcamp = rtcamp;
                vector<double> uc_init(mx * my, 0);
                uc_init = uc;
                // vector<int> pulse_marker_c(N, 0);

                // find camp source cell location
                // always let the frontmost red cell be the source Oct19 2022
                double src_cell[2] = {0, 0};
                double max_x = 0;
                for (int i = 0; i < N; i++)
                {
                    if (X[i][0] >= max_x && LB[i] == 0)
                    {
                        src_cell[0] = X[i][0];
                        src_cell[1] = X[i][1];
                    }
                }

                // find the grid point closest to the source
                vector<int> src_loc_c{0, 0};
                double dist = Lx;
                for (int j = 0; j < mx; j++)
                {
                    for (int k = 0; k < my; k++)
                    {
                        double gridP[2] = {x[j], y[k]};
                        double tmp_dist = 0;
                        for (int l = 0; l < 2; l++)
                        {
                            tmp_dist += (src_cell[l] - gridP[l]) * (src_cell[l] - gridP[l]);
                        }
                        tmp_dist = sqrt(tmp_dist);
                        if (tmp_dist < dist)
                        {
                            dist = tmp_dist;
                            src_loc_c[0] = j;
                            src_loc_c[1] = k;
                        }
                    }
                }
                FCpara.src_loc = src_loc_c;
                // solve camp nonlinear pde
                sr_init_c = sr_out_c;
                ADI_nonlin(uc_init, sr_init_c, tstart, tend, X, LB, tcamp, PCpara, FCpara, uc, rtcamp, sr_out_c, pulse_marker_c);

                /*****************************************************************
                                    Cell Differentiation
                            compute cumulated differentiation feedback
                ******************************************************************/
                for (int j = 0; j < mx; j++)
                {
                    for (int k = 0; k < my; k++)
                    {
                        int pos = k * mx + j;
                        rud[j][k] = ud[pos];
                        ruc[j][k] = uc[pos];
                    }
                }
                // start collecting feedback info
                if (time_count >= Tlb[lb_check] - int_fdb_t && time_count <= Tlb[lb_check])
                {
                    // concentrations at cell locations
                    vector<double> ruc_x(N, 0);
                    vector<double> rud_x(N, 0);
                    for (int i = 0; i < N; i++)
                    {
                        double cellP[2] = {0, 0};
                        cellP[0] = Xmod[i][0];
                        cellP[1] = Xmod[i][1]; // the i-th cell's location
                        double dist = Lx;
                        int k_clo = 0;
                        int j_clo = 0;
                        for (int k = 0; k < my; k++)
                        {
                            for (int j = 0; j < mx; j++)
                            {
                                double gridP[2] = {x[j], y[k]};
                                double tmp_dist = 0;
                                for (int l = 0; l < 2; l++)
                                {
                                    tmp_dist += (cellP[l] - gridP[l]) * (cellP[l] - gridP[l]);
                                }
                                tmp_dist = sqrt(tmp_dist);
                                if (tmp_dist < dist)
                                {
                                    dist = tmp_dist;
                                    k_clo = k;
                                    j_clo = j;
                                }
                            }
                        }
                        ruc_x[i] = ruc[j_clo][k_clo];
                        rud_x[i] = rud[j_clo][k_clo];
                    }
                    // integrate the feedback function
                    for (int i = 0; i < N; i++)
                    {
                        Feedback_pt[i] += (pow(rud_x[i], exp_d) / (pow(a0, exp_d) + pow(rud_x[i], exp_d))) * (pow(b1, exp_c) / (pow(b1, exp_c) + pow(ruc_x[i], exp_c))) * timestep / int_fdb_t;
                        Feedback_tp[i] += (pow(b0, exp_d) / (pow(b0, exp_d) + pow(rud_x[i], exp_d))) * (pow(ruc_x[i], exp_c) / (pow(a1, exp_c) + pow(ruc_x[i], exp_c))) * timestep / int_fdb_t;
                    }
                }

            } // end of time loop
            // cout << "solve cell motion equations successfully" << "\n";

            /*****************************************************************
                                    Update Variables from
                                motion and signalling equations
            ******************************************************************/
            // write output into files
            // only extract the last generation of output
            for (int i = 0; i < (N + Tnc) * 2; i++)
            {
                if (i % 2 == 0)
                {
                    *foutm << NV_Ith_S(init_cond, i) << " " << NV_Ith_S(init_cond, i + 1) << "\n";
                }
            }
            for (int j = 0; j < mx; j++)
            {
                for (int k = 0; k < my; k++)
                {
                    *foutd << rud[j][k] << " ";
                    *foutc << ruc[j][k] << " ";
                }
                *foutd << "\n";
                *foutc << "\n";
            }

            for (int i = 0; i < N; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    int startP = N + pre_count(i, j, nc, N);
                    int endP = startP + nc[i][j];
                    int count = 0;
                    for (int k = startP; k < endP; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            Sc[i][j][count][l] = NV_Ith_S(init_cond, 2 * k + l); // init_cond[k][l]
                        }
                        count += 1;
                    }
                }
            }

            /*****************************************************************
                                        Chemotaxis
                            (update direction of camp argC)
            ******************************************************************/

            vector<double> ux(N, 0);
            vector<double> uy(N, 0); // partial derivative
            for (int i = 0; i < N; i++)
            {
                double cellP[2] = {0, 0};
                cellP[0] = Xmod[i][0];
                cellP[1] = Xmod[i][1]; // the i-th cell's location
                double dist = Lx;
                int k_clo = 0;
                int j_clo = 0;
                for (int k = 0; k < my; k++)
                {
                    for (int j = 0; j < mx; j++)
                    {
                        double gridP[2] = {x[j], y[k]};
                        double tmp_dist = 0;
                        for (int l = 0; l < 2; l++)
                        {
                            tmp_dist += (cellP[l] - gridP[l]) * (cellP[l] - gridP[l]);
                        }
                        tmp_dist = sqrt(tmp_dist);
                        if (tmp_dist < dist)
                        {
                            dist = tmp_dist;
                            k_clo = k;
                            j_clo = j;
                        }
                    }
                } //  find the one closest to cell i
                if (j_clo == 0)
                {
                    ux[i] = (ruc[j_clo + 1][k_clo] - ruc[j_clo][k_clo]) / dx;
                }
                else if (j_clo == Nx)
                {
                    ux[i] = (ruc[j_clo][k_clo] - ruc[j_clo - 1][k_clo]) / dx;
                }
                if (k_clo == 0)
                {
                    uy[i] = (ruc[j_clo][k_clo + 1] - ruc[j_clo][k_clo]) / dy;
                }
                else if (k_clo == Ny)
                {
                    uy[i] = (ruc[j_clo][k_clo] - ruc[j_clo][k_clo - 1]) / dy;
                }
                if (j_clo != 0 && k_clo != 0 && j_clo != Nx && k_clo != Ny)
                {
                    ux[i] = (ruc[j_clo + 1][k_clo] - ruc[j_clo - 1][k_clo]) / (2 * dx);
                    uy[i] = (ruc[j_clo][k_clo + 1] - ruc[j_clo][k_clo - 1]) / (2 * dy);
                }
                ArgC[i] = atan2(uy[i], ux[i]); // lie in [-pi,pi]
            }
            gc_para.argC = ArgC;

            /*****************************************************************
                                Cell Differentiation
                                (update label variable LB)
            ******************************************************************/

            if (time_count >= Tlb[lb_check])
            { // check cell type every eight minutes
                lb_check += 1;
                cout << "check cell type at t = " << time_count << "\n";

                // check the cumulated feedback
                for (int i = 0; i < N; i++)
                {
                    if (LB[i] == 1 && Feedback_pt[i] >= tsh_PT && Feedback_tp[i] < tsh_TP)
                    {
                        triggerState[i] = 1;
                        vector<int> diff_weights;
                        diff_weights.push_back(1);
                        diff_weights.push_back(19);
                        std::discrete_distribution<> diff(diff_weights.begin(), diff_weights.end());
                        LB[i] = diff(generator);
                        if (LB[i] == 0)
                        {
                            cout << "cell " << i << " turn into PST"
                                 << "\n";
                        }
                    }
                    else if (LB[i] == 0 && Feedback_pt[i] < tsh_PT && Feedback_tp[i] >= tsh_TP)
                    {
                        triggerState[i] = 1;
                        vector<int> diff_weights;
                        diff_weights.push_back(9);
                        diff_weights.push_back(1);
                        std::discrete_distribution<> diff(diff_weights.begin(), diff_weights.end());
                        LB[i] = diff(generator);
                        if (LB[i] == 1)
                        {
                            cout << "cell " << i << " turn into PSP"
                                 << "\n";
                        }
                    }
                    else if (Feedback_pt[i] >= tsh_PT && Feedback_tp[i] >= tsh_TP)
                    {
                        triggerState[i] = 2;
                        cout << "cell " << i << " turn into PST or PSP"
                             << "\n";
                        vector<int> diff_weights;
                        diff_weights.push_back(1);
                        diff_weights.push_back(1);
                        std::discrete_distribution<> diff(diff_weights.begin(), diff_weights.end());
                        LB[i] = diff(generator);
                    }
                    else
                    {
                        triggerState[i] = 0;
                        vector<int> diff_weights;
                        diff_weights.push_back(19); // tend to stay as the same
                        diff_weights.push_back(1);
                        std::discrete_distribution<> diff(diff_weights.begin(), diff_weights.end());
                        int change_or_not = diff(generator);
                        if (change_or_not == 1)
                        {
                            LB[i] = 1 - LB[i];
                        }
                    }
                }
                // clear all cumulated feedback
                std::fill(Feedback_pt.begin(), Feedback_pt.end(), 0);
                std::fill(Feedback_tp.begin(), Feedback_tp.end(), 0);
            } // check cell type

            for (int i = 0; i < N; i++)
            {
                *foutm << triggerState[i] << " ";
            }
            *foutm << "\n";
        } // if time_count >= 20

        /*****************************************************************
                               Execute Actions
                             (update Sc, nc, Tc)
        ******************************************************************/
        // remove out-of-range cadherins and integrins --------------------
        for (int i = 0; i < N; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                for (int k = 0; k < nc[i][j]; k++)
                {
                    double d1 = 0;
                    double d2 = 0;
                    for (int l = 0; l < 2; l++)
                    {
                        d1 += (X[i][l] - Sc[i][j][k][l]) * (X[i][l] - Sc[i][j][k][l]);
                        d2 += (X[j][l] - Sc[i][j][k][l]) * (X[j][l] - Sc[i][j][k][l]);
                    }
                    d1 = sqrt(d1);
                    d2 = sqrt(d2);
                    if (d1 > r + Lc || d2 > r + Lc)
                    {
                        int tbr_idx[3] = {i, j, k};
                        cad.removeC(Sc, nc, Tc, tbr_idx, N, Nc);
                    }
                }
            }
        }

        // execute "action"s rc ac ----------------------------------------
        if (action == "rc")
        {
            // find the cadherin to be removed
            int tbr_i = -1;
            int tbr_j = -1;
            int tbr_k = -1;
            for (int i = 0; i < N; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    for (int k = 0; k < nc[i][j]; k++)
                    {
                        if (Tc[i][j][k] == t1)
                        {
                            tbr_i = i;
                            tbr_j = j;
                            tbr_k = k;
                        }
                    }
                }
            }
            // remove the cadherin
            if (tbr_i != -1 && tbr_j != -1 && tbr_k != -1)
            { // if it was not already removed
                int tbr_idx[3] = {tbr_i, tbr_j, tbr_k};
                cout << "Remove the c-site " << Sc[tbr_i][tbr_j][tbr_k][0] << " , " << Sc[tbr_i][tbr_j][tbr_k][1]
                     << " between cell " << tbr_i << " and cell " << tbr_k << "\n";
                cad.removeC(Sc, nc, Tc, tbr_idx, N, Nc);
            }
            // advance Tc by t0
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < nc[i][j]; k++)
                    {
                        Tc[i][j][k] += -t0;
                    }
                }
            }
        }
        else if (action == "ac")
        {
            // Tc should -t0
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < nc[i][j]; k++)
                    {
                        Tc[i][j][k] += -t0;
                    }
                }
            }
            vector<vector<int>> nnc(N, vector<int>(N, 0));                                                                                      // int n_c[N * N];
            vector<vector<vector<vector<double>>>> nSc(N, vector<vector<vector<double>>>(N, vector<vector<double>>(Nc, vector<double>(2, 0)))); // double S_c[Tnc][2];
            vector<vector<vector<double>>> nTc(N, vector<vector<double>>(N, vector<double>(Nc, 0)));                                            // double T_c[Tnc];
            cad.gen_c(X, gc_para, nSc, nnc, nTc, LB);

            vector<double> pullout = {};
            int count_nonzero = 0;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < Nc; k++)
                    {
                        if (nTc[i][j][k] > 0)
                        {
                            double nonzero_select = nTc[i][j][k];
                            pullout.push_back(nonzero_select);
                            count_nonzero += 1;
                        }
                    }
                }
            }

            if (count_nonzero != 0)
            { // if there exists some new cadherins
                // randomly select one c-site
                double Psample = 0;
                Psample = *select_randomly(pullout.begin(), pullout.end());

                // find the c-site
                int tba_i = 0;
                int tba_j = 0;
                int tba_k = 0;
                for (int i = 0; i < N; i++)
                {
                    for (int j = i + 1; j < N; j++)
                    {
                        for (int k = 0; k < Nc; k++)
                        {
                            if (nTc[i][j][k] == Psample)
                            {
                                tba_i = i;
                                tba_j = j;
                                tba_k = k;
                                cout << "the new c-site: " << nSc[tba_i][tba_j][tba_k][0] << " "
                                     << nSc[tba_i][tba_j][tba_k][1] << "\n";
                            }
                        }
                    }
                }

                int cidx_add[2] = {tba_i, tba_j};                                                // add between cell i and cell j
                double ccor_add[2] = {nSc[tba_i][tba_j][tba_k][0], nSc[tba_i][tba_j][tba_k][1]}; // coordiates of the new c-site
                double ctime_add = nTc[tba_i][tba_j][tba_k];                                     // duration of the new c-site
                cout << "Add a new c-site " << ccor_add[0] << " , " << ccor_add[1]
                     << " between cell " << tba_i << " and cell " << tba_j << "\n";
                cad.addC(Sc, nc, Tc, cidx_add, ccor_add, ctime_add, N, Nc);
            }
        }

        /*****************************************************************
                            Restore attachment sites
                              (avoid falling apart)
        ******************************************************************/
        if (time_count >= 1)
        { // after the initial noise goes away
            for (int i = 0; i < N; i++)
            {
                // check total # of c-sites for cell i
                int Tnc_single = 0;
                for (int j = 0; j < N; j++)
                {
                    Tnc_single += nc[i][j]; // sum up the i-th row in nc
                }

                if (Tnc_single < mNc)
                { // if the total # of c-sites of cell i is below mNc
                    // all the candidates within the cell i reach-out range
                    vector<int> cand_j = {};
                    int cand_count = 0;
                    for (int j = 0; j < N; j++)
                    {
                        double dis = 0;
                        for (int l = 0; l < 2; l++)
                        {
                            dis += (X[i][l] - X[j][l]) * (X[i][l] - X[j][l]);
                        }
                        dis = sqrt(dis);
                        if (dis > 0 && dis <= 2 * r + 2 * Lc)
                        {
                            cand_j.push_back(j);
                            cand_count += 1;
                        }
                    }

                    if (cand_count != 0)
                    { // if there exists available neighbors
                        for (int k = 0; k < mNc - Tnc_single; k++)
                        { // try (mNc - Tnc_single) times to form c-sites with neighboring cells
                            // randomly select a candidate
                            int j = 0;
                            j = *select_randomly(cand_j.begin(), cand_j.end());

                            // form a new c-site
                            vector<int> prob_weights;
                            int prob = -1;
                            if (LB[i] == 0 && LB[j] == 0)
                            { // TT
                                prob_weights.push_back(1);
                                prob_weights.push_back(19);
                                std::discrete_distribution<> d(prob_weights.begin(), prob_weights.end());
                                prob = d(generator);
                            }
                            else if (LB[i] == 1 && LB[j] == 1)
                            { // PP
                                prob_weights.push_back(1);
                                prob_weights.push_back(19);
                                std::discrete_distribution<> d(prob_weights.begin(), prob_weights.end());
                                prob = d(generator);
                            }
                            else
                            { // TP
                                prob_weights.push_back(1);
                                prob_weights.push_back(19);
                                std::discrete_distribution<> d(prob_weights.begin(), prob_weights.end());
                                prob = d(generator);
                            }
                            if (prob == 0)
                            { // half the time will add the c-site back
                                int cidx_add[2] = {i, j};
                                double ccor_add[2] = {0, 0};
                                double ctime_add = 0;
                                cad.gen_c_single(X, i, j, gc_para, ccor_add, ctime_add, LB);
                                if (ccor_add[0] != 0 && ccor_add[1] != 0)
                                {
                                    cout << "Restore a new c-site " << ccor_add[0] << " , " << ccor_add[1]
                                         << " between cell " << i << " and cell " << j << "\n";
                                    cad.addC(Sc, nc, Tc, cidx_add, ccor_add, ctime_add, N, Nc);
                                }
                            }
                        }
                    }
                }
            }
        } // restore attachment sites

        time_count += t0;
        cout << "time count: " << time_count << "\n";

        // deallocate solution vector
        N_VDestroy(init_cond);
        // clear ode memmory since the next generation has different size
        CVodeFree(&cvode_mem);
    } // end of one generation (outermost while loop)

    // free memory  ---------------------------------------------------------------
    SUNContext_Free(&sunctx);
    foutm->close();
    foutc->close();
    foutd->close();
    delete foutm;
    delete foutc;
    delete foutd;

    getchar();
    return 0;
} // end of main

/************************************************************************************
                    auxillary function definitions
************************************************************************************/
double norm(double *x, int dim)
{
    // calculate the euclidean norm of a vector
    double norm = 0;
    for (int i = 0; i < dim; i++)
    {
        norm += x[i] * x[i];
    }
    norm = sqrt(norm);
    return norm;
}

double sigmoid(double x, double m, double c, double s)
{
    double y = 0;
    y = m / (1 + exp(-(x - c) / s));
    return y;
}

void redim21(double (*V)[2], int len_V, double *v)
{
    // reshape the two dimensional array V into one dimensional array v
    // len_V is the length of V (to be converted)
    for (int i = 0; i < len_V * 2; i++)
    {
        if (i % 2 == 0)
        {
            v[i] = V[i / 2][0];
            v[i + 1] = V[i / 2][1];
        }
    }
}

// -------------------------------------------------------------------------------
void redim12(double *v, int len_v, double (*V)[2])
{
    // reshape the one dimensional array v into two dimensional array V
    // len_v is the length of v
    for (int i = 0; i < len_v; i++)
    {
        if (i % 2 == 0)
        {
            V[i / 2][0] = v[i];
            V[i / 2][1] = v[i + 1];
        }
    }
}

// -------------------------------------------------------------------------------
int pre_count(int i, int j, vector<vector<int>> &nc, int N)
{
    // add all the number in the upper half of nc
    // before the (i,j) entry of nc (i,j start from 0)
    // N is the size of the matrix nc

    // only consider the part above the diagonal (i < j)
    if (i > j)
    {
        int tmp = i;
        i = j;
        j = tmp;
    }

    int y = 0;
    if (i == 0)
    { // the first row of nc
        for (int jj = 0; jj <= j - 1; jj++)
        {
            y += nc[i][jj];
        }
    }
    else
    {
        for (int ii = 0; ii <= i - 1; ii++)
        {
            for (int jj = ii + 1; jj <= N - 1; jj++)
            {
                y += nc[ii][jj];
            }
        }
        for (int jj = i + 1; jj <= j - 1; jj++)
        {
            y += nc[i][jj];
        }
    }
    return y;
}

// -------------------------------------------------------------------------------
void Id_c(int nij, double *xi, double *xj, vector<vector<double>> &sij, double r, double L, vector<double> &Phi)
{
    for (int i = 0; i < nij; i++)
    {
        Phi[i] = 0;
    }
    for (int k = 0; k < nij; k++)
    {
        double d1 = 0;
        double d2 = 0;
        double norms = 0;
        for (int l = 0; l < 2; l++)
        {
            d1 += (xi[l] - sij[k][l]) * (xi[l] - sij[k][l]);
            d2 += (xj[l] - sij[k][l]) * (xj[l] - sij[k][l]);
            norms += sij[k][l] * sij[k][l];
        }
        d1 = sqrt(d1);
        d2 = sqrt(d2);
        if (d1 <= r + L && d2 <= r + L && norms != 0)
        {
            Phi[k] = 1;
        }
    }
}

// -------------------------------------------------------------------------------
void f_blue(int nij, double alpha, double r, double lc, double Lc, double *x, vector<vector<double>> &s, vector<vector<double>> &F)
{
    for (int i = 0; i < nij; i++)
    {
        for (int l = 0; l < 2; l++)
        {
            F[i][l] = 0;
        }
    }
    for (int k = 0; k < nij; k++)
    {
        double d = 0;
        for (int l = 0; l < 2; l++)
        {
            d += (x[l] - s[k][l]) * (x[l] - s[k][l]);
        }
        d = sqrt(d);
        if (d > r + lc && d < r + Lc)
        { // only drag forces, no pushaway forces
            for (int l = 0; l < 2; l++)
            {
                F[k][l] = alpha * (d - r - lc) * (x[l] - s[k][l]) / d;
            }
        }
    }
}

// -------------------------------------------------------------------------------
void f_red(int nij, double alpha, double r, double lc, double Lc, double nl, double *x, vector<vector<double>> &s, vector<vector<double>> &F)
{
    for (int i = 0; i < nij; i++)
    {
        for (int l = 0; l < 2; l++)
        {
            F[i][l] = 0;
        }
    }
    for (int k = 0; k < nij; k++)
    {
        double d = 0;
        for (int l = 0; l < 2; l++)
        {
            d += (x[l] - s[k][l]) * (x[l] - s[k][l]);
        }
        d = sqrt(d);
        if (x[0] >= s[k][0])
        { // if the cadherin is formed in the back
            if (d < r + Lc)
            {
                for (int l = 0; l < 2; l++)
                {
                    F[k][l] = alpha * (d - nl * (r + lc)) * (x[l] - s[k][l]) / d;
                }
            }
        }
        else if (x[0] < s[k][0])
        { // if the cadherin is formed in the front
            if (d > r + lc && d < r + Lc)
            { // only drag forces, no pushaway forces
                for (int l = 0; l < 2; l++)
                {
                    F[k][l] = alpha * (d - r - lc) * (x[l] - s[k][l]) / d;
                }
            }
        }
    }
}

// -------------------------------------------------------------------------------
void f1_blue(double alpha, double r, double lc, double Lc, double *x, double *y, double *F)
{
    // force exerted by blue cell on the cadherin
    // the nij = 1 version of f_blue
    for (int l = 0; l < 2; l++)
    {
        F[l] = 0;
    }
    double d = 0;
    for (int l = 0; l < 2; l++)
    {
        d += (x[l] - y[l]) * (x[l] - y[l]);
    }
    d = sqrt(d);
    if (d > r + lc && d < r + Lc)
    {
        for (int l = 0; l < 2; l++)
        {
            F[l] = alpha * (d - r - lc) * (x[l] - y[l]) / d;
        }
    }
}

// -------------------------------------------------------------------------------
void f1_red(double alpha, double r, double lc, double Lc, double nl, double *x, double *y, double *F)
{
    // force exerted by red cell on the cadherin
    // the nij = 1 version of f_red
    for (int l = 0; l < 2; l++)
    {
        F[l] = 0;
    }
    double d = 0;
    for (int l = 0; l < 2; l++)
    {
        d += (x[l] - y[l]) * (x[l] - y[l]);
    }
    d = sqrt(d);
    if (x[0] <= y[0])
    { // when the cadherin is ahead of the cell
        if (d > r + lc && d < r + Lc)
        {
            for (int l = 0; l < 2; l++)
            {
                F[l] = alpha * (d - r - lc) * (x[l] - y[l]) / d;
            }
        }
    }
    else if (x[0] > y[0])
    { // when the cadherin is behind the cell
        if (d < r + Lc)
        {
            for (int l = 0; l < 2; l++)
            {
                F[l] = alpha * (d - nl * (r + lc)) * (x[l] - y[l]) / d;
            }
        }
    }
}

// -------------------------------------------------------------------------------
double Id_b(double *xi, double *xj, double r, double lc)
{
    double Phi = 0;

    double d = 0;
    for (int l = 0; l < 2; l++)
    {
        d += (xi[l] - xj[l]) * (xi[l] - xj[l]);
    }
    d = sqrt(d);

    if (d > 0 && d < 2 * r + 2 * lc)
    {
        Phi = 1;
    }

    return Phi;
}

// -------------------------------------------------------------------------------
void g(double alpha, double alphaM, double r, double lc, double *x, double *y, double *B)
{
    for (int l = 0; l < 2; l++)
    {
        B[l] = 0;
    }

    double d = 0;
    for (int l = 0; l < 2; l++)
    {
        d += (x[l] - y[l]) * (x[l] - y[l]);
    }
    d = sqrt(d); // norm of x-y

    if (d > 2 * r && d < 2 * r + 2 * lc)
    {
        for (int l = 0; l < 2; l++)
        {
            B[l] = alphaM * (2 * r + 2 * lc - d) * (x[l] - y[l]) / d;
        }
    }
    else if (d > 0 && d < 2 * r)
    {
        for (int l = 0; l < 2; l++)
        {
            B[l] = alphaM * (2 * lc) * (x[l] - y[l]) / d +
                   alpha * (exp(2 * (2 * r - d)) - 1) * (x[l] - y[l]) / d;
        }
    }
}

// -------------------------------------------------------------------------------
int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        fprintf(stderr,
                "\n SUNDIALS_ERROR : %s()  failed  -  returned   NULL   pointer \n\n",
                funcname);
        return (1);
    }
    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *)flagvalue;
        if (*errflag < 0)
        {
            fprintf(stderr,
                    "\n SUNDIALS_ERROR : %s()  failed   with   flag  = %d \n\n",
                    funcname, *errflag);
            return (1);
        }
    }
    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        fprintf(stderr,
                "\n MEMORY_ERROR : %s() failed  -  returned   NULL pointer \n\n", funcname);
        return (1);
    }
    return (0);
}

// -------------------------------------------------------------------------------
int f_rhs(double t, N_Vector y, N_Vector ydot, void *user_data)
{
    // typedef int (*CVRhsFn)(realtype t, N Vector y, N Vector ydot, void *user_data);
    //  This function computes the ODE right-hand side for a given value of the independent variable t and state vector y.
    //  t is the current value of the independent variable.
    //  y is the current value of the dependent variable vector, y(t).
    //  ydot is the output vector f(t; y).
    //  user_data is the user data pointer passed to CVodeSetUserData.
    //  A CVRhsFn should return 0 if successful, a positive value if a recoverable error occurred (in which case cvode will attempt to correct),
    //  or a negative value if it failed unrecoverably (in which case the integration is halted and CV RHSFUNC FAIL is returned).
    //  Allocation of memory for ydot is handled within cvode.

    // user data
    Mpara *mpara;
    mpara = (Mpara *)user_data;

    // extract parameters
    int N = mpara->N; // int Ni = mpara->Ni;
    double r = mpara->r;
    double lc = mpara->lc;
    double Lc = mpara->Lc;
    double nl = mpara->nl;
    double alphaB = mpara->alphaB;
    double alphaS = mpara->alphaS;
    vector<double> AlphaM(3, 0);
    AlphaM = mpara->AlphaM;
    double mu = mpara->mu;
    double nu = mpara->nu;
    vector<double> AlphaC(2, 0);
    AlphaC = mpara->AlphaC;
    vector<int> LB(N, 0);
    LB = mpara->label;
    vector<vector<int>> nc(N, vector<int>(N, 0));
    nc = mpara->nc;

    // copy current value y(t) N_Vector --> V(t) array
    int Tnc = pre_count(N - 1, N, nc, N); // total number of c-sites

    vector<double> V((N + Tnc) * 2, 0); // double V[(N + Pnc) * 2];
    for (int i = 0; i < (N + Tnc) * 2; i++)
    {
        V[i] = NV_Ith_S(y, i);
    }

    // reshape current value V(t) --> rV(t)
    double V_array[10000];
    for (int i = 0; i < (N + Tnc) * 2; i++)
    {
        V_array[i] = V[i];
    }
    double rV[5000][2];
    redim12(V_array, (N + Tnc) * 2, rV);

    // Mc encodes forces due to cell-cell connection (c-sites)
    vector<vector<vector<double>>> Mc(N, vector<vector<double>>(N, vector<double>(2, 0))); // double Mc[N][N][2];

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // check cell type and determine alphaC
            double alphaC = 0;
            if (LB[i] == 0)
            {
                alphaC = AlphaC[0];
            }
            else if (LB[i] == 1)
            {
                alphaC = AlphaC[1];
            }

            int nc_ij = nc[i][j]; // number of c-sites between i and j
            double Vi[2];
            double Vj[2];
            for (int l = 0; l < 2; l++)
            {
                Vi[l] = rV[i][l];
                Vj[l] = rV[j][l];
            }

            if (nc_ij != 0)
            { // if there exists c-sites between i and j
                // grab all the c-sites between i and j
                // Create a vector containing n vectors of size m.
                vector<vector<double>> Sc(nc_ij, vector<double>(2, 0)); // double Sc[nc_ij][2];
                int startP = N + pre_count(i, j, nc, N);
                int endP = startP + nc_ij;
                int count = 0;
                for (int k = startP; k < endP; k++)
                {
                    for (int l = 0; l < 2; l++)
                    {
                        Sc[count][l] = rV[k][l];
                    }
                    count += 1;
                }
                // the forces due to c-sites
                vector<vector<double>> Fc(nc_ij, vector<double>(2, 0));
                // asymmetric membrane for red cells
                if (LB[i] == 0 && LB[j] == 1)
                {
                    f_red(nc_ij, alphaC, r, lc, Lc, nl, Vi, Sc, Fc);
                }
                else
                {
                    f_blue(nc_ij, alphaC, r, lc, Lc, Vi, Sc, Fc);
                }
                for (int l = 0; l < 2; l++)
                {
                    for (int k = 0; k < nc_ij; k++)
                    {
                        Mc[i][j][l] += Fc[k][l];
                    }
                }
            }
        }
    }

    // G encodes body forces between cells
    vector<vector<vector<double>>> G(N, vector<vector<double>>(N, vector<double>(2, 0))); // double G[N][N][2];

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double Vi[2];
            double Vj[2];
            for (int l = 0; l < 2; l++)
            {
                Vi[l] = rV[i][l];
                Vj[l] = rV[j][l];
            }
            if (j != i)
            {
                double alphaM = 0;
                if (LB[i] == 0 && LB[j] == 0)
                { // PST PST
                    alphaM = AlphaM[0];
                }
                else if (LB[i] == 1 && LB[j] == 1)
                { // PSP PSP
                    alphaM = AlphaM[1];
                }
                else
                { // PST PSP
                    alphaM = AlphaM[2];
                }
                double Phig_c = 0;
                double Bc[2] = {0, 0};
                Phig_c = Id_b(Vi, Vj, r, lc);
                g(alphaB, alphaM, r, lc, Vi, Vj, Bc);
                for (int l = 0; l < 2; l++)
                {
                    G[i][j][l] = Bc[l] * Phig_c;
                }
            }
        }
    }

    // B encodes body forces from the substrate
    vector<vector<double>> B(N, vector<double>(2, 0));
    for (int i = 0; i < N; i++)
    {
        double foot[2] = {0, 0};
        foot[0] = rV[i][0];
        double Vi[2];
        for (int l = 0; l < 2; l++)
        {
            Vi[l] = rV[i][l];
        }
        double Phig_o = 0;
        double Bo[2] = {0, 0};
        Phig_o = Id_b(Vi, foot, r, li);
        double alphaM = 0;
        g(alphaS, alphaM, r, lc, Vi, foot, Bo);
        for (int l = 0; l < 2; l++)
        {
            B[i][l] = Bo[l] * Phig_o;
        }
    }

    // W encodes body forces from the upper wall
    vector<vector<double>> W(N, vector<double>(2, 0));
    double wall_y = 8;
    for (int i = 0; i < N; i++)
    {
        double foot[2] = {0, wall_y};
        foot[0] = rV[i][0];
        double Vi[2];
        for (int l = 0; l < 2; l++)
        {
            Vi[l] = rV[i][l];
        }
        double Phig_o = 0;
        double Wo[2] = {0, 0};
        Phig_o = Id_b(Vi, foot, r, 0);
        double alphaM = 0;
        g(alphaS, alphaM, r, 0, Vi, foot, Wo);
        for (int l = 0; l < 2; l++)
        {
            W[i][l] = Wo[l] * Phig_o;
        }
    }

    // ode system
    vector<vector<double>> dVdt(N + Tnc, vector<double>(2, 0));

    for (int i = 0; i < N; i++)
    {
        double sum_Mc[2] = {0, 0};
        double sum_G[2] = {0, 0};
        for (int l = 0; l < 2; l++)
        {
            for (int j = 0; j < N; j++)
            {
                sum_Mc[l] += Mc[i][j][l];
                sum_G[l] += G[i][j][l];
            }
        }
        for (int l = 0; l < 2; l++)
        {
            dVdt[i][l] = -sum_Mc[l] / mu + sum_G[l] / mu + B[i][l] / mu + W[i][l] / mu;
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            // check cell type and determine alphaC
            double alphaC_1 = 0;
            double alphaC_2 = 0;
            if (LB[i] == 0 && LB[j] == 0)
            {
                alphaC_1 = AlphaC[0];
                alphaC_2 = AlphaC[0];
            }
            else if (LB[i] == 1 && LB[j] == 1)
            {
                alphaC_1 = AlphaC[1];
                alphaC_2 = AlphaC[1];
            }
            else if (LB[i] == 0 && LB[j] == 1)
            {
                alphaC_1 = AlphaC[0];
                alphaC_2 = AlphaC[1];
            }
            else if (LB[i] == 1 && LB[j] == 0)
            {
                alphaC_1 = AlphaC[1];
                alphaC_2 = AlphaC[0];
            }

            double Vi[2];
            double Vj[2];
            for (int l = 0; l < 2; l++)
            {
                Vi[l] = rV[i][l];
                Vj[l] = rV[j][l];
            }
            int startP = N + pre_count(i, j, nc, N);
            int endP = startP + nc[i][j];
            for (int k = startP; k < endP; k++)
            {
                double Vk[2];
                for (int l = 0; l < 2; l++)
                {
                    Vk[l] = rV[k][l]; // c-sites between i and j
                }
                double Fi[2];
                double Fj[2];
                if (LB[i] == 0 && LB[j] == 1)
                {
                    f1_red(alphaC_1, r, lc, Lc, nl, Vi, Vk, Fi);
                    f1_blue(alphaC_2, r, lc, Lc, Vj, Vk, Fj);
                }
                else if (LB[i] == 1 && LB[j] == 0)
                {
                    f1_blue(alphaC_1, r, lc, Lc, Vi, Vk, Fi);
                    f1_red(alphaC_2, r, lc, Lc, nl, Vj, Vk, Fj);
                }
                else
                {
                    f1_blue(alphaC_1, r, lc, Lc, Vi, Vk, Fi);
                    f1_blue(alphaC_2, r, lc, Lc, Vj, Vk, Fj);
                }
                // motion equation for this cadherin
                for (int l = 0; l < 2; l++)
                {
                    dVdt[k][l] = Fi[l] / nu + Fj[l] / nu;
                }
            }
        }
    }

    // reshape dVdt --> rdVdt
    double dVdt_array[5000][2];
    for (int i = 0; i < 5000; i++)
    {
        for (int l = 0; l < 2; l++)
        {
            dVdt_array[i][l] = 0;
        }
    }
    for (int i = 0; i < N + Tnc; i++)
    {
        for (int l = 0; l < 2; l++)
        {
            dVdt_array[i][l] = dVdt[i][l];
        }
    }
    double rdVdt[10000];
    redim21(dVdt_array, N + Tnc, rdVdt);

    // rdVdt --> ydot output
    // ydot = N_VMake_Serial((N + Tnc) * 2, rdVdt); // creates a new N_Vector instead of modify ydot
    // NV_LENGTH_S(ydot) = (N + Tnc) * 2;
    for (int i = 0; i < (N + Tnc) * 2; i++)
    {
        NV_Ith_S(ydot, i) = rdVdt[i];
    }
    return (0);
}

// -------------------------------------------------------------------------------
void gen_init(int N, vector<vector<double>> &X, vector<vector<int>> &nc, vector<vector<vector<vector<double>>>> &Sc, double *cc_array)
{
    // this function generate a vector cc of initial conditions
    // based on the current cell locations X and c-sites locations Sc

    // total number of c-sites
    int Tnc = pre_count(N - 1, N, nc, N);

    // pack locations of cell centers and c-sites
    vector<vector<double>> Tmpx(N + Tnc, vector<double>(2, 0)); // double Tmpx[N + Tnc][2];
    for (int i = 0; i < N; i++)
    {
        for (int l = 0; l < 2; l++)
        {
            Tmpx[i][l] = X[i][l];
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            int startP = N + pre_count(i, j, nc, N);
            int endP = N + pre_count(i, j, nc, N) + nc[i][j];
            int count = 0;
            for (int k = startP; k < endP; k++)
            {
                for (int l = 0; l < 2; l++)
                {
                    Tmpx[k][l] = Sc[i][j][count][l]; // some of these c-sites are zero
                }
                count += 1; // in the end count should = nc[i][j]
            }
        }
    }

    // reshape Tmpx into one dimension array tmpx
    double Tmpx_array[5000][2];
    if (N + Tnc > 5000)
    {
        cout << "error: N + Tnc > 5000"
             << "\n";
    }
    for (int i = 0; i < N + Tnc; i++)
    {
        for (int l = 0; l < 2; l++)
        {
            Tmpx_array[i][l] = Tmpx[i][l];
        }
    }

    // set vector of initial values
    redim21(Tmpx_array, N + Tnc, cc_array);
}

// -------------------------------------------------------------------------------
void ADI_nonlin(vector<double> &u_init, vector<double> &sr_init, double tstart, double tend, vector<vector<double>> &X, vector<int> label, vector<double> &tchem, PDEpara CKpara, NLpara Fpara, vector<double> &u, vector<double> &rtchem, vector<double> &sr_out, vector<int> &pulse_marker)
{
    // u_init: the initial conditions
    // tstart: initial time
    // tend: advance solution
    // X: cell locations
    // label: cell types of each cell
    // tchem: the time xi output the chemical last time
    // CKpara: parameters for CK algorithm
    // Fpara: parameters for nonlinear forcing term

    int N = CKpara.N;
    double Lx = CKpara.Lx;
    double Ly = CKpara.Ly;
    int Nx = CKpara.Nx;
    int Ny = CKpara.Ny;
    double mindt = CKpara.mindt;
    double sigma = CKpara.sigma;
    double tol = CKpara.tol;

    double baseT = CKpara.baseT;
    double baseP = CKpara.baseP;
    double tshT = CKpara.tshT;
    double tshP = CKpara.tshP;
    double Tt = CKpara.Tt;
    double Tp = CKpara.Tp;
    double slopT = CKpara.slopT;
    double slopP = CKpara.slopP;

    double T = tend - tstart;
    double dx = Lx / Nx;
    double dy = Ly / Ny;

    // determine time step
    int Nt = 0;
    double dt = 0;
    if (T < mindt)
    {
        Nt = 1;
        dt = T;
    }
    else
    {
        Nt = ceil(T / mindt);
        dt = T / Nt;
    }
    // cout << "number of time steps: " << Nt << "\n";

    // number of points including boundray points
    int mx = Nx + 1;
    int my = Ny + 1;

    // grid points including boundaries
    vector<double> x(mx, 0);
    vector<double> y(my, 0);
    for (int j = 0; j < mx; j++)
    {
        x[j] = j * dx;
    }
    for (int k = 0; k < my; k++)
    {
        y[k] = k * dy;
    }

    // parameters
    double rx = (sigma * dt) / (dx * dx);
    double thetax = sigma / (dx * dx);
    double ry = (sigma * dt) / (dy * dy);
    double thetay = sigma / (dy * dy);

    // matrices ---------------------------------------------
    // diagonal vector of Cx
    for (int i = 0; i < Nx - 1; i++)
    {
        Dx[i] = 1 + rx;
    }
    // subdiagonal vector of Cx
    for (int i = 0; i < Nx - 2; i++)
    {
        Bx[i] = -rx / 2;
    }
    // diagonal vector of Cy
    for (int j = 0; j < Ny - 1; j++)
    {
        Dy[j] = 1 + ry;
    }
    // subdiagonal vector of Cy
    for (int j = 0; j < Ny - 2; j++)
    {
        By[j] = -ry / 2;
    }

    // initialization ----------------------------------------
    u = u_init;
    vector<double> uh(mx * my, 0);
    uh = u;
    vector<double> upre(mx * my, 0);
    vector<double> F(mx * my, 0);
    vector<double> Fpre(mx * my, 0);

    // starting stage of sr
    vector<double> sr(N, 0);
    sr = sr_init;
    nonlin_F(u_init, sr, X, label, Fpara, x, y, F);

    // iteratively solve the linear system -------------------------------------------------
    for (int n = 0; n < Nt; n++)
    {
        // vectors
        upre = u; // u_{n-1}
        Fpre = F; // F(u_{n-1})

        // matrices
        vector<vector<double>> U(mx, vector<double>(my, 0));
        vector<vector<double>> FF(mx, vector<double>(my, 0));
        for (int j = 0; j < mx; j++)
        {
            for (int k = 0; k < my; k++)
            {
                int pos = k * mx + j;
                U[j][k] = u[pos];
                FF[j][k] = F[pos];
            }
        }

        /************************************************************************************
                                        ADI method
        ************************************************************************************/
        int INFO = 0;
        int NRHS = 1;
        int orderAx = Nx - 1;
        int LDBx = Nx - 1;
        int orderAy = Ny - 1;
        int LDBy = Ny - 1;

        // predictor step: from u_{n-1} obtain u_{n-1/2} -------------------------------------------------------
        vector<vector<double>> UH(mx, vector<double>(my, 0)); // initialize u_{n-1/2}
        vector<double> rhs1(Nx - 1, 0);
        float RHS1[NX - 1] = {0}; // array version of rhs1
        for (int k = 1; k < Ny; k++)
        {
            for (int j = 1; j < Nx; j++)
            {
                RHS1[j - 1] = ry * U[j][k - 1] / 2 + (1 - ry) * U[j][k] + ry * U[j][k + 1] / 2 + FF[j][k];
            }
            sptsv_(&orderAx, &NRHS, Dx, Bx, RHS1, &LDBx, &INFO);
            // solve the system A * X = B where A is symmetric positive definite tridiagonal
            // X and B are N-by-NRHS matrices
            // N: integer, the order of matrix A
            // NRHS: the number of right hand sides
            // D: real array, dimension N. On entry, the n diagonal elements of A; On exit, the n diagonal elements of D from A = L*D*L**T
            // E: real aaray, dimension N-1. On entry, the n-1 subdiagonal elements of A; On exit, the n-1 subdiagonal elements of L
            // B: real array, dimension (LDB, NRHS). On entry, the N-by-nrhs right hand side matrix. On exit, if INFO = 0, the solution matrix X
            // LDB: integer, leading dimension of array B
            // INFO: integer, =0 successful exit
            for (int j = 1; j < Nx; j++)
            {
                UH[j][k] = RHS1[j - 1];
                // if (UH[j][k] == 0){
                //     cout << "concentration UH is zero" << "\n";
                // }
                // else{
                //     cout << "concentration UH is nonzero" << "\n";
                // }
                // if (j == 14 && k == 10 && UH[j][k] == U[j][k]){
                //     cout << "UH = Upre at cell 1" << "\n";
                // }
            }
        }
        // turn matrix UH into vector uh
        for (int j = 0; j < mx; j++)
        {
            for (int k = 0; k < my; k++)
            {
                int pos = k * mx + j;
                uh[pos] = UH[j][k];
            }
        }

        // corrector step: from u_{n-1/2} obtain u_n -------------------------------------------
        vector<double> rhs2(Ny - 1, 0);
        float RHS2[NY - 1] = {0}; // array version of rhs1
        for (int j = 1; j < Nx; j++)
        {
            for (int k = 1; k < Ny; k++)
            {
                RHS2[k - 1] = rx * UH[j - 1][k] / 2 + (1 - rx) * UH[j][k] + rx * UH[j + 1][k] / 2 + FF[j][k];
            }
            sptsv_(&orderAy, &NRHS, Dy, By, RHS2, &LDBy, &INFO);
            for (int k = 1; k < Ny; k++)
            {
                U[j][k] = RHS2[k - 1];
            }
        }
        // turn matrix U into vector u
        for (int j = 0; j < mx; j++)
        {
            for (int k = 0; k < my; k++)
            {
                int pos = k * mx + j;
                u[pos] = (U[j][k] > 0) ? U[j][k] : 0;
            }
        }

        // nonlinear term at u_{n} ----------------------------------------------------------------------
        // construct previous ut at t_{n-1/2}
        vector<double> ut_pre(mx * my, 0);
        for (int j = 0; j < mx; j++)
        {
            for (int k = 0; k < my; k++)
            {
                // ut_pre[i] = (u[i] - uh[i]) / (dt / 2);
                int pos = k * mx + j;
                ut_pre[pos] = (u[pos] - upre[pos]) / dt;
                // if (u[pos] == upre[pos] && j == 14 && k == 10){
                //     cout << "concentration does not change at cell 1" << "\n";
                // }

                // ut_pre does not change means u - upre does not change
            }
        }
        // reshape ut_pre as a grid function
        vector<vector<double>> rut_pre(mx, vector<double>(my, 0));
        for (int j = 0; j < mx; j++)
        {
            for (int k = 0; k < my; k++)
            {
                int pos = k * mx + j;
                rut_pre[j][k] = ut_pre[pos];
            }
        }

        // find values of ut at locations of xi
        vector<vector<double>> Xmod(N, vector<double>(2, 0));
        for (int i = 0; i < N; i++)
        {
            Xmod[i][0] = std::fmod(X[i][0], Lx);
            Xmod[i][1] = std::fmod(X[i][1], Ly);
        }
        vector<double> ut_pre_x(N, 0);
        for (int i = 0; i < N; i++)
        {
            double dist = Lx; // changed location on Dec 8 !!!!!
            double cellP[2] = {0, 0};
            cellP[0] = Xmod[i][0];
            cellP[1] = Xmod[i][1]; // the i-th cell's location
            int k_clo = 0;
            int j_clo = 0;
            for (int k = 0; k < my; k++)
            {
                for (int j = 0; j < mx; j++)
                {
                    double gridP[2] = {x[j], y[k]};
                    double tmp_dist = 0;
                    for (int l = 0; l < 2; l++)
                    {
                        tmp_dist += (cellP[l] - gridP[l]) * (cellP[l] - gridP[l]);
                    }
                    tmp_dist = sqrt(tmp_dist);

                    if (tmp_dist < dist)
                    {
                        dist = tmp_dist;
                        k_clo = k;
                        j_clo = j;
                        // cout << "intermediate: " << j_clo << ", " << k_clo << " is closer to cell " << i << "\n";
                    }
                }
            }
            if (k_clo == 0 && j_clo == 0)
            {
                cout << "the approximate grid for cell " << i << " is (0,0)"
                     << "\n";
            }
            ut_pre_x[i] = rut_pre[j_clo][k_clo];

            // cout stimulation information for each cell, inside the Nt time loop
            // if (label[i] == 0){
            //     if (ut_pre_x[i] > tshT){
            //         cout << "cell " << i << " is stimulated" << "\n";
            //     }
            //     else{
            //         cout << "ut_pre_x[" << i << "] = " << ut_pre_x[i] << "\n";
            //     }
            // }
            // if (label[i] == 1){
            //     if (ut_pre_x[i] > tshP){
            //         cout << "cell " << i << " is stimulated" << "\n";
            //     }
            //     else{
            //         cout << "ut_pre_x[" << i << "] = " << ut_pre_x[i] << "\n";
            //     }
            // }
            // cout information end ---------------------------------------------
        }

        // construct chemical-output function at location xi
        vector<double> tdiff(N, 0);
        for (int i = 0; i < N; i++)
        {
            tdiff[i] = tstart + n * dt - tchem[i];
        } // time past since xi output chemical last time

        for (int i = 0; i < N; i++)
        {
            if (label[i] == 0)
            {
                if (ut_pre_x[i] > tshT && tdiff[i] >= Tt)
                {
                    tchem[i] = tstart + n * dt;
                    tdiff[i] = tstart + n * dt - tchem[i];
                }
                if (tdiff[i] <= Tt / 4)
                {
                    sr[i] = baseT + slopT * tdiff[i];
                }
                else if (tdiff[i] > Tt / 4 && tdiff[i] <= Tt / 2)
                {
                    sr[i] = slopT * Tt / 2 + baseT - slopT * tdiff[i];
                }
                else
                {
                    sr[i] = baseT;
                }
            }
            else if (label[i] == 1)
            {
                if (ut_pre_x[i] > tshP && tdiff[i] >= Tp)
                {
                    tchem[i] = tstart + n * dt;
                    tdiff[i] = tstart + n * dt - tchem[i];
                }

                if (tdiff[i] <= Tp / 4)
                {
                    sr[i] = baseP + slopP * tdiff[i];
                }
                else if (tdiff[i] > Tp / 4 && tdiff[i] <= Tp / 2)
                {
                    sr[i] = slopP * Tp / 2 + baseP - slopP * tdiff[i];
                }
                else
                {
                    sr[i] = baseP;
                }
            }
        }

        // construct nonlinear term
        nonlin_F(u, sr, X, label, Fpara, x, y, F);

        // mark pulsing cells in the last Nt time loop
        if (n == Nt - 1)
        {
            for (int i = 0; i < N; i++)
            {
                if (label[i] == 0)
                {
                    if (ut_pre_x[i] > tshT)
                    {
                        // cout << "cell " << i << " is stimulated" << "\n";
                        pulse_marker[i] = 1; // stimulated PST cells
                    }
                    // else{
                    //     cout << "ut_pre_x[" << i << "] = " << ut_pre_x[i] << "\n";
                    // }
                }
                if (label[i] == 1)
                {
                    if (ut_pre_x[i] > tshP)
                    {
                        // cout << "cell " << i << " is stimulated" << "\n";
                        pulse_marker[i] = 1; // stimulated PSP cells
                    }
                    // else{
                    //     cout << "ut_pre_x[" << i << "] = " << ut_pre_x[i] << "\n";
                    // }
                }
            }
        }
    }
    rtchem = tchem;
    sr_out = sr;
}

// -------------------------------------------------------------------------------
void nonlin_F(vector<double> &upre, vector<double> &sr, vector<vector<double>> &X, vector<int> &label, NLpara para, vector<double> &x, vector<double> &y, vector<double> &F)
{
    // upre: previous u at gridpoints
    // ut_pre: previous ut at cell center locations
    // tc: the time when xi output camp last time

    int N = para.N;
    int Nx = para.Nx;
    int Ny = para.Ny;
    int K = para.K;
    double Lx = para.Lx;
    double Ly = para.Ly;
    double src = para.src;
    vector<int> src_loc{0, 0};
    src_loc = para.src_loc;
    double Vc = para.Vc;
    double V0 = para.V0;
    double gm6_t = para.gm6_t;
    double gm6_p = para.gm6_p;
    double gm7_t = para.gm7_t;
    double gm7_p = para.gm7_p;
    double gm8 = para.gm8;
    double gm9 = para.gm9;

    double dx = Lx / Nx;
    double dy = Ly / Ny;
    double mx = Nx + 1;
    double my = Ny + 1;
    double mu = 0;
    double sigma = dx;
    if (dy > dx)
    {
        sigma = dy;
    }

    // degrade everywhere
    for (int i = 0; i < mx * my; i++)
    {
        F[i] = -gm9 * upre[i] / (upre[i] + gm8);
    }

    // mode cell lacations
    vector<vector<double>> Xmod(N, vector<double>(2, 0));
    for (int i = 0; i < N; i++)
    {
        Xmod[i][0] = std::fmod(X[i][0], Lx);
        Xmod[i][1] = std::fmod(X[i][1], Ly);
    }

    for (int k = 0; k < my; k++)
    {
        for (int j = 0; j < mx; j++)
        {
            double gridP[2] = {x[j], y[k]};
            for (int i = 0; i < N; i++)
            {
                double dist = 0;
                for (int l = 0; l < 2; l++)
                {
                    dist += (X[i][l] - gridP[l]) * (X[i][l] - gridP[l]);
                }
                dist = sqrt(dist);
                double delta = exp(-(dist - mu) * (dist - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
                int position = k * mx + j;
                double coeff = K * Vc * delta / V0;
                if (label[i] == 0)
                {
                    F[position] += coeff * (sr[i] - gm7_t * upre[position] / (upre[position] + gm6_t));
                }
                else if (label[i] == 1)
                {
                    F[position] += coeff * (sr[i] - gm7_p * upre[position] / (upre[position] + gm6_p));
                }
            }
        }
    }

    int Src_Loc = src_loc[1] * mx + src_loc[0];
    F[Src_Loc] = F[Src_Loc] + src;
}
