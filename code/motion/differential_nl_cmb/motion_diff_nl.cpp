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
                                  constants
************************************************************************************/
double tol = 1e-4;                    // tolerance for steady states
double T_MAX = 70;                    // maximum time for motion equations
double uni_dt = 1e-5;                 // uniform delta t
double r = 1;                         // radius of the nucleus (micrometer)
double lc = 2;                        // natural length of the spring for c-sites
double li = 1;                        // natural length of the spring for i-sites
double Lc = 8;                        // max length of the spring for c-sites (micrometer)
double nl = 2;                        // diff natural length
double Li = 7;                        // max length of the spring for i-sites (micrometer)
double argC = 0;                      // angle of the gradient of cAMP concentration [-pi,pi]
double rangeC = pi / 16;              // reach out range
int N = 80;                           // total number of cells
int rowN = 20;                        // number of cells in each row
int columnN = 4;                      // number of cells in each column
int Nc_max = 4;                       // max number of c-sites formed between two cells
vector<int> Nc(N, Nc_max);            // number of attempts to form cadherins for each cell
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

int main()
{
    // read in nl values
    std::ifstream fin;
    fin.open("nl_value.dat");
    fin >> nl;

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
    gc_para.Nc_max = Nc_max;
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

    // time partitions -----------------------------------------
    // form c/i-sites frequently
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

    // create output files ------------------------------------------
    string filename1 = "motion.dat";
    ofstream *foutm = new ofstream(filename1);

    // initialize cell centers
    vector<vector<double>> X(N, vector<double>(2, 0));
    double gap = 2 * r;
    // initial rectangular shape
    for (int i = 0; i < N; i++)
    {
        int posx = fmod(i, rowN);
        int posy = i / rowN;
        X[i][0] = r + lc + posx * gap;
        X[i][1] = r + lc + posy * gap;
    }

    // initialize cell types
    vector<int> LB(N, 1); // 0: PST; 1: PSP
    LB[20] = 0;
    LB[60] = 0;
    LB[1] = 0;
    LB[21] = 0;
    LB[41] = 0;
    LB[22] = 0;
    LB[42] = 0;
    LB[62] = 0;
    LB[3] = 0;
    LB[23] = 0;
    LB[63] = 0;
    LB[24] = 0;
    LB[44] = 0;
    LB[5] = 0;
    LB[25] = 0;
    LB[45] = 0;

    // initialize all the c-sites ------------------------------
    vector<double> ArgC(N, argC);
    vector<vector<int>> nc(N, vector<int>(N, 0));                                                                                          // int n_c[N * N];
    vector<vector<vector<vector<double>>>> Sc(N, vector<vector<vector<double>>>(N, vector<vector<double>>(Nc_max, vector<double>(2, 0)))); // double S_c[Tnc][2];
    vector<vector<vector<double>>> Tc(N, vector<vector<double>>(N, vector<double>(Nc_max, 0)));                                            // double T_c[Tnc];
    gc_para.argC = ArgC;

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
                                    motion system
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

        // the c-site that wears out the fastest
        double t1 = 100;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < Nc_max; k++)
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

        // write number of cadherins for this generation
        *foutm << Tnc << "\n";

        // write how long this fixed group of c-sites & i-sites last
        *foutm << t0 << "\n";

        // wirte cell labels
        for (int i = 0; i < N; i++)
        {
            *foutm << LB[i] << " ";
        }
        *foutm << "\n";

        // run cell motion equations, advance t0 ----------------------------------------------
        if (time_count >= 20)
        {
            for (int ti = 0; ti < 100; ti++)
            {
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

                // write output cells and c-sites locations
                cout << "time: " << tstart << "\n";
                for (int i = 0; i < (N + Tnc) * 2; i++)
                {
                    if (i % 2 == 0)
                    {
                        *foutm << NV_Ith_S(init_cond, i) << " " << NV_Ith_S(init_cond, i + 1) << "\n";
                    }
                }
            } // end of time loop
            cout << "solve cell motion equations successfully"
                 << "\n";

            // update cell centers X --------------------------------------
            for (int i = 0; i < 2 * N; i++)
            {
                if (i % 2 == 0)
                {
                    X[i / 2][0] = NV_Ith_S(init_cond, i);
                    X[i / 2][1] = NV_Ith_S(init_cond, i + 1);
                }
            }

            // update c-sites locations Sc --------------------------------
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
        }

        // update Sc, nc, Tc -----------------------------------------------------
        // remove out-of-range c-sites
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
                        cad.removeC(Sc, nc, Tc, tbr_idx, N, Nc_max);
                    }
                }
            }
        }

        // execute "action"s rc ac -----------------------------------------------------
        if (action == "rc")
        {
            // find the c-sites to be removed
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

            // remove the c-site
            if (tbr_i != -1 && tbr_j != -1 && tbr_k != -1)
            { // if it was not already removed
                int tbr_idx[3] = {tbr_i, tbr_j, tbr_k};
                cout << "Remove the c-site " << Sc[tbr_i][tbr_j][tbr_k][0] << " , " << Sc[tbr_i][tbr_j][tbr_k][1]
                     << " between cell " << tbr_i << " and cell " << tbr_k << "\n";
                cad.removeC(Sc, nc, Tc, tbr_idx, N, Nc_max);
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
            vector<vector<int>> nnc(N, vector<int>(N, 0));                                                                                          // int n_c[N * N];
            vector<vector<vector<vector<double>>>> nSc(N, vector<vector<vector<double>>>(N, vector<vector<double>>(Nc_max, vector<double>(2, 0)))); // double S_c[Tnc][2];
            vector<vector<vector<double>>> nTc(N, vector<vector<double>>(N, vector<double>(Nc_max, 0)));                                            // double T_c[Tnc];
            cad.gen_c(X, gc_para, nSc, nnc, nTc, LB);

            vector<double> pullout = {};
            int count_nonzero = 0;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < Nc_max; k++)
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
                        for (int k = 0; k < Nc_max; k++)
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
                cad.addC(Sc, nc, Tc, cidx_add, ccor_add, ctime_add, N, Nc_max);
            }
        }

        // restore c-sites and i-sites -------------------------------------------------------
        if (time_count >= 1)
        { // after the initial noise goes away
            // check minimum number of c-sites -----------------------------
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
                                    cad.addC(Sc, nc, Tc, cidx_add, ccor_add, ctime_add, N, Nc_max);
                                }
                            }
                        }
                    }
                }
            }
        }

        // complete all updates and move on -------------------------------------------
        time_count += t0;
        cout << "time count: " << time_count << "\n";

        // deallocate solution vector
        N_VDestroy(init_cond);
        // clear ode memmory since the next generation has different size
        CVodeFree(&cvode_mem);
    } // end of one generation for a fixed group of c-sites

    // free memory  ---------------------------------------------------------------
    SUNContext_Free(&sunctx);
    foutm->close();
    delete foutm;

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
    double wall_y = 10;
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
