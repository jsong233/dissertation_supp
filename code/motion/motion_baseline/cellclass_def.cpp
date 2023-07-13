#include "cellclass.h"
#include <cmath>
#include <cfloat>
#include <fstream>
#include <random>

/************************************************************************************
                                  Cadherin
************************************************************************************/
void Cadherin::gen_c(vector<vector<double>> &X, GCpara gc_para, vector<vector<vector<vector<double>>>> &Sc, vector<vector<int>> &nc, vector<vector<vector<double>>> &Tc, vector<int> &LB)
{
    // randomly generate new c-sites from cells X
    std::random_device rd;
    std::mt19937 generator(rd());

    int N = gc_para.N;
    int Nc_max = gc_para.Nc_max;
    vector<int> Nc(N, 0);
    Nc = gc_para.Nc;
    double r = gc_para.r;
    double Lc = gc_para.Lc;
    vector<double> argC(N, 0);
    argC = gc_para.argC;
    vector<double> rangeC(N, 0);
    rangeC = gc_para.rangeC;
    double tc_PP = gc_para.tc_PP;
    double tc_TT = gc_para.tc_TT;
    double tc_TP = gc_para.tc_TP;

    vector<vector<int>> nnc(N, vector<int>(N, 0));                                                                                         // int nnc[N][N];            // intermediate number of c-sites
    vector<vector<vector<vector<double>>>> sc(N, vector<vector<vector<double>>>(N, vector<vector<double>>(Nc_max, vector<double>(2, 0)))); // double sc[N][N][Nc][2];      // intermediate possible c-sites
    // double Sc[N][N][Nc][2];      // all the possible c-sites
    // double Tc[N][N][Nc];         // attach time for each c-site
    // int nc[N][N];                // pairwise number of c-sites

    // generate all possible c-sites in a pairwise manner
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double d = 0;
            for (int l = 0; l < 2; l++)
            {
                d += (X[i][l] - X[j][l]) * (X[i][l] - X[j][l]);
            }
            d = sqrt(d); // distance between cell i and cell j

            if (j != i && d <= 2 * r + 2 * Lc)
            {
                // the vector pointing from i to j
                double delta[2];
                for (int l = 0; l < 2; l++)
                {
                    delta[l] = X[j][l] - X[i][l];
                }

                for (int k = 0; k < Nc[i]; k++)
                { // Nc[i] attempts for cell i

                    // outreach length unformly distributed between r and Lc
                    double br = 0;
                    std::uniform_real_distribution<> Udistr(r, r + Lc);
                    br = Udistr(generator);

                    // outreach angle normally distributed
                    std::normal_distribution<> Ndistr(argC[i], rangeC[i]);
                    double btheta = Ndistr(generator);

                    // the outreach vector b and the candidate cadherin
                    double b[2] = {br * cos(btheta), br * sin(btheta)};
                    double sc_cand[2] = {X[i][0] + b[0], X[i][1] + b[1]};

                    // whether the outreach vector b is within a reasonable range of delta
                    // whether the outreach direction aligns with the cell ij direction
                    // find cos(angle) between b and delta
                    double dot_b_delta = 0;
                    double norm_b = 0;
                    double norm_delta = 0;
                    for (int l = 0; l < 2; l++)
                    {
                        dot_b_delta += b[l] * delta[l];
                        norm_b += b[l] * b[l];
                        norm_delta += delta[l] * delta[l];
                    }
                    norm_b = sqrt(norm_b);
                    norm_delta = sqrt(norm_delta);

                    // whether the candidate cadherin is also close enough to cell j
                    double dis_to_j = 0;
                    for (int l = 0; l < 2; l++)
                    {
                        dis_to_j += (X[j][l] - sc_cand[l]) * (X[j][l] - sc_cand[l]);
                    }
                    dis_to_j = sqrt(dis_to_j);

                    // combine criterion 1 and criterion 2
                    // if the angle between b and delta <= pi/8
                    // and the candidate cadherin is within the outreach range of cell j
                    if (dot_b_delta / (norm_b * norm_delta) >= cos(pi / 8) && dis_to_j <= r + Lc)
                    {
                        nnc[i][j] += 1;               // add one c-site between i,j
                        int position = nnc[i][j] - 1; // the k-th c-site locate at k-1
                        for (int l = 0; l < 2; l++)
                        {
                            sc[i][j][position][l] = sc_cand[l];
                        }
                    }
                } // end of for k
            }     // end of if
        }
    }

    // truncate c-sites
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            int count = nnc[i][j] + nnc[j][i];
            if (count > Nc_max)
            {
                if (nnc[j][i] < nnc[i][j])
                {
                    nnc[j][i] = Nc_max - nnc[i][j];
                }
                else
                {
                    nnc[i][j] = Nc_max - nnc[j][i];
                }
                //   nnc[i][j] = Nc * nnc[i][j] / count;
                //   nnc[j][i] = Nc - nnc[i][j];
            } // make sure nnc_ij + nnc_ji <= Nc
        }
    }

    // concatenate c-sites
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            int count = nnc[i][j] + nnc[j][i];
            for (int k = 0; k < count; k++)
            {
                if (k < nnc[i][j])
                {
                    for (int l = 0; l < 2; l++)
                    {
                        Sc[i][j][k][l] = sc[i][j][k][l];
                    }
                }
                else
                {
                    for (int l = 0; l < 2; l++)
                    {
                        int posLeft = k - nnc[i][j];
                        Sc[i][j][k][l] = sc[j][i][posLeft][l];
                    }
                }
            }
            nc[i][j] = nnc[i][j] + nnc[j][i];
            nc[j][i] = nc[i][j];
            for (int k = 0; k < Nc_max; k++)
            {
                for (int l = 0; l < 2; l++)
                {
                    Sc[j][i][k][l] = Sc[i][j][k][l];
                }
            }
        }
    }

    // attach time
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            for (int k = 0; k < nc[i][j]; k++)
            {
                double time = 0;
                while (time <= 0)
                {
                    if (LB[i] == 0 && LB[j] == 0)
                    {
                        std::normal_distribution<> Ndistr(tc_TT, 0.05); // deviation 0.2 * mean
                        time = Ndistr(generator);
                    }
                    else if (LB[i] == 1 && LB[j] == 1)
                    {
                        std::normal_distribution<> Ndistr(tc_PP, 0.05); // deviation 0.2 * mean
                        time = Ndistr(generator);
                    }
                    else
                    {
                        std::normal_distribution<> Ndistr(tc_TP, 0.05); // deviation 0.2 * mean
                        time = Ndistr(generator);
                    }
                }
                Tc[i][j][k] = time;
                Tc[j][i][k] = Tc[i][j][k];
            }
        }
    }

    // Sc[N][N][Nc][2];
    // Tc[N][N][Nc];
    // nc[N][N];
}

// --------------------------------------------------------------------------
void Cadherin::gen_c_init(vector<vector<double>> &X, GCpara gc_para, vector<vector<vector<vector<double>>>> &Sc, vector<vector<int>> &nc, vector<vector<vector<double>>> &Tc, vector<int> &LB)
{
    // randomly generate new c-sites from cells X
    std::random_device rd;
    std::mt19937 generator(rd());

    int N = gc_para.N;
    int Nc_max = gc_para.Nc_max;
    Nc_max = Nc_max / 2;
    vector<int> Nc(N, 0);
    Nc = gc_para.Nc;
    // initialize with fewer c-sites
    for (int i = 0; i < N; i++)
    {
        Nc[i] = Nc[i] / 2;
    }
    double r = gc_para.r;
    double Lc = gc_para.Lc;
    vector<double> argC(N, 0);
    argC = gc_para.argC;
    vector<double> rangeC(N, 0);
    rangeC = gc_para.rangeC;
    double tc_PP = gc_para.tc_PP;
    double tc_TT = gc_para.tc_TT;
    double tc_TP = gc_para.tc_TP;

    vector<vector<int>> nnc(N, vector<int>(N, 0));                                                                                         // int nnc[N][N];            // intermediate number of c-sites
    vector<vector<vector<vector<double>>>> sc(N, vector<vector<vector<double>>>(N, vector<vector<double>>(Nc_max, vector<double>(2, 0)))); // double sc[N][N][Nc][2];      // intermediate possible c-sites
                                                                                                                                           // double Sc[N][N][Nc][2];      // all the possible c-sites
                                                                                                                                           // double Tc[N][N][Nc];         // attach time for each c-site
                                                                                                                                           // int nc[N][N];                // pairwise number of c-sites

    // generate all possible c-sites in a pairwise manner
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double d = 0;
            for (int l = 0; l < 2; l++)
            {
                d += (X[i][l] - X[j][l]) * (X[i][l] - X[j][l]);
            }
            d = sqrt(d); // distance between cell i and cell j

            if (j != i && d <= 2 * r + 2 * Lc)
            {
                // the vector pointing from i to j
                double delta[2];
                for (int l = 0; l < 2; l++)
                {
                    delta[l] = X[j][l] - X[i][l];
                }

                for (int k = 0; k < Nc[i]; k++)
                {

                    // outreach length unformly distributed between r and Lc
                    double br = 0;
                    std::uniform_real_distribution<> Udistr(r, r + Lc);
                    br = Udistr(generator); // unformly distributed between r and Lc

                    // outreach angle normally distributed
                    std::normal_distribution<> Ndistr(argC[i], rangeC[i]);
                    double btheta = Ndistr(generator);

                    // the outreach vector b and the candidate cadherin
                    double b[2] = {br * cos(btheta), br * sin(btheta)};
                    double sc_cand[2] = {X[i][0] + b[0], X[i][1] + b[1]};

                    // whether the outreach vector b is within a reasonable range of delta
                    // whether the outreach direction aligns with the cell ij direction
                    // find cos(angle) between b and delta
                    double dot_b_delta = 0;
                    double norm_b = 0;
                    double norm_delta = 0;
                    for (int l = 0; l < 2; l++)
                    {
                        dot_b_delta += b[l] * delta[l];
                        norm_b += b[l] * b[l];
                        norm_delta += delta[l] * delta[l];
                    }
                    norm_b = sqrt(norm_b);
                    norm_delta = sqrt(norm_delta);

                    // whether the candidate cadherin is also close enough to cell j
                    double dis_to_j = 0;
                    for (int l = 0; l < 2; l++)
                    {
                        dis_to_j += (X[j][l] - sc_cand[l]) * (X[j][l] - sc_cand[l]);
                    }
                    dis_to_j = sqrt(dis_to_j);

                    // combine criterion 1 and criterion 2
                    // if the angle between b and delta <= pi/8
                    // and the candidate cadherin is within the outreach range of cell j
                    if (dot_b_delta / (norm_b * norm_delta) >= cos(pi / 8) && dis_to_j <= r + Lc)
                    {
                        nnc[i][j] += 1;               // add one c-site between i,j
                        int position = nnc[i][j] - 1; // the k-th c-site locate at k-1
                        for (int l = 0; l < 2; l++)
                        {
                            sc[i][j][position][l] = sc_cand[l];
                        }
                    }
                } // end of for k
            }     // end of if
        }
    }

    // truncate c-sites
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            int count = nnc[i][j] + nnc[j][i];
            if (count > Nc_max)
            {
                if (nnc[j][i] < nnc[i][j])
                {
                    nnc[j][i] = Nc_max - nnc[i][j];
                }
                else
                {
                    nnc[i][j] = Nc_max - nnc[j][i];
                }
            } // make sure nnc_ij + nnc_ji <= Nc
        }
    }

    // concatenate c-sites
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            int count = nnc[i][j] + nnc[j][i];
            for (int k = 0; k < count; k++)
            {
                if (k < nnc[i][j])
                {
                    for (int l = 0; l < 2; l++)
                    {
                        Sc[i][j][k][l] = sc[i][j][k][l];
                    }
                }
                else
                {
                    for (int l = 0; l < 2; l++)
                    {
                        int posLeft = k - nnc[i][j];
                        Sc[i][j][k][l] = sc[j][i][posLeft][l];
                    }
                }
            }
            nc[i][j] = nnc[i][j] + nnc[j][i];
            nc[j][i] = nc[i][j];
            for (int k = 0; k < Nc_max; k++)
            {
                for (int l = 0; l < 2; l++)
                {
                    Sc[j][i][k][l] = Sc[i][j][k][l];
                }
            }
        }
    }

    // attach time
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            for (int k = 0; k < nc[i][j]; k++)
            {
                double time = 0;
                while (time <= 0)
                {
                    if (LB[i] == 0 && LB[j] == 0)
                    {
                        std::uniform_real_distribution<> Udistr(tc_TT / 2, tc_TT * 2);
                        time = Udistr(generator);
                    }
                    else if (LB[i] == 1 && LB[j] == 1)
                    {
                        std::uniform_real_distribution<> Udistr(tc_PP / 2, tc_PP * 2);
                        time = Udistr(generator);
                    }
                    else
                    {
                        std::uniform_real_distribution<> Udistr(tc_TP / 2, tc_TP * 2);
                        time = Udistr(generator);
                    }
                }
                Tc[i][j][k] = time;
                Tc[j][i][k] = Tc[i][j][k];
            }
        }
    }

    // Sc[N][N][Nc][2];
    // Tc[N][N][Nc];
    // nc[N][N];
}

// --------------------------------------------------------------------------
void Cadherin::gen_c_single(vector<vector<double>> &X, int tba_i, int tba_j, GCpara gc_para, double *Sc, double &Tc, vector<int> &LB)
{
    // generate a single c-site between cell i and cell j
    // Sc: coordinates of the new c-site
    // Tc: attach time of the new c-site

    std::random_device rd;
    std::mt19937 generator(rd());

    int N = gc_para.N;
    int Nc_max = gc_para.Nc_max;
    vector<int> Nc(N, 0);
    Nc = gc_para.Nc;
    double r = gc_para.r;
    double Lc = gc_para.Lc;
    vector<double> argC(N, 0);
    argC = gc_para.argC;
    vector<double> rangeC(N, 0);
    rangeC = gc_para.rangeC;
    double tc_PP = gc_para.tc_PP;
    double tc_TT = gc_para.tc_TT;
    double tc_TP = gc_para.tc_TP;
    int i = tba_i;
    int j = tba_j;
    vector<vector<double>> sc(2, vector<double>(2, 0)); // save two possible c-sites

    // distance between cell i and cell j
    double d = 0;
    for (int l = 0; l < 2; l++)
    {
        d += (X[i][l] - X[j][l]) * (X[i][l] - X[j][l]);
    }
    d = sqrt(d);

    // generate a c-site between celli and cellj
    if (d <= 2 * r + 2 * Lc)
    {
        // the vector pointing from i to j
        double delta[2];
        for (int l = 0; l < 2; l++)
        {
            delta[l] = X[j][l] - X[i][l];
        }

        // outreach length unformly distributed between r and Lc
        std::uniform_real_distribution<> Udistr(r, Lc);
        double br = Udistr(generator);

        // outreach angle normally distributed
        std::normal_distribution<> Ndistr(argC[i], rangeC[i]);
        double btheta = Ndistr(generator);

        // the outreach vector b and the candidate cadherin
        double b[2] = {br * cos(btheta), br * sin(btheta)};
        double sc_cand[2] = {X[i][0] + b[0], X[i][1] + b[1]};

        // whether the outreach vector b is within a reasonable range of delta
        // whether the outreach direction aligns with the cell ij direction
        // find cos(angle) between b and delta
        double dot_b_delta = 0;
        double norm_b = 0;
        double norm_delta = 0;
        for (int l = 0; l < 2; l++)
        {
            dot_b_delta += b[l] * delta[l];
            norm_b += b[l] * b[l];
            norm_delta += delta[l] * delta[l];
        }
        norm_b = sqrt(norm_b);
        norm_delta = sqrt(norm_delta);

        // whether the candidate cadherin is also close enough to cell j
        double dis_to_j = 0;
        for (int l = 0; l < 2; l++)
        {
            dis_to_j += (X[j][l] - sc_cand[l]) * (X[j][l] - sc_cand[l]);
        }
        dis_to_j = sqrt(dis_to_j);

        // combine criterion 1 and criterion 2
        // if the angle between b and delta <= pi/8
        // and the candidate cadherin is within the outreach range of cell j
        if (dot_b_delta / (norm_b * norm_delta) >= cos(pi / 8) && dis_to_j <= r + Lc)
        { // change to cos(pi/8) on Jul29
            for (int l = 0; l < 2; l++)
            {
                sc[0][l] = sc_cand[l];
            }
        }
    }

    // generate a c-site between cellj and celli
    if (d <= 2 * r + 2 * Lc)
    {
        // the vector pointing from j to i
        double delta[2];
        for (int l = 0; l < 2; l++)
        {
            delta[l] = X[j][l] - X[i][l];
        }

        // outreach length unformly distributed between r and Lc
        std::uniform_real_distribution<> Udistr(r, r + Lc);
        double br = Udistr(generator); // unformly distributed between r and Lc

        // outreach angle normally distributed
        std::normal_distribution<> Ndistr(argC[j], rangeC[i]);
        double btheta = Ndistr(generator);

        // the outreach vector b and the candidate cadherin
        double b[2] = {br * cos(btheta), br * sin(btheta)};
        double sc_cand[2] = {X[j][0] + b[0], X[j][1] + b[1]};

        // whether the outreach vector b is within a reasonable range of delta
        // whether the outreach direction aligns with the cell ij direction
        // find cos(angle) between b and delta
        double dot_b_delta = 0;
        double norm_b = 0;
        double norm_delta = 0;
        for (int l = 0; l < 2; l++)
        {
            dot_b_delta += b[l] * delta[l];
            norm_b += b[l] * b[l];
            norm_delta += delta[l] * delta[l];
        }
        norm_b = sqrt(norm_b);
        norm_delta = sqrt(norm_delta);

        // whether the candidate cadherin is also close enough to cell i
        double dis_to_i = 0;
        for (int l = 0; l < 2; l++)
        {
            dis_to_i += (X[i][l] - sc_cand[l]) * (X[i][l] - sc_cand[l]);
        }
        dis_to_i = sqrt(dis_to_i);

        // combine criterion 1 and criterion 2
        // if the angle between b and delta <= pi/8
        // and the candidate cadherin is within the outreach range of cell i
        if (dot_b_delta / (norm_b * norm_delta) >= cos(pi / 8) && dis_to_i <= r + Lc)
        {
            for (int l = 0; l < 2; l++)
            {
                sc[1][l] = sc_cand[l];
            }
        }
    }

    // choose betweewn two possible c-sites
    double check0 = sc[0][0] * sc[0][0] + sc[0][1] * sc[0][1];
    double check1 = sc[1][0] * sc[1][0] + sc[1][1] * sc[1][1];
    if (check0 == 0 && check1 != 0)
    {
        Sc[0] = sc[1][0];
        Sc[1] = sc[1][1];
    }
    else if (check0 != 0 && check1 == 0)
    {
        Sc[0] = sc[0][0];
        Sc[1] = sc[0][1];
    }
    else if (check0 != 0 && check1 != 0)
    {
        vector<int> weights = {1, 1};
        std::discrete_distribution<> d(weights.begin(), weights.end());
        int target = d(generator);
        Sc[0] = sc[target][0];
        Sc[1] = sc[target][1];
    }

    // attach time
    double time = 0;
    while (time <= 0)
    {
        if (LB[tba_i] == 0 && LB[tba_j] == 0)
        {
            std::normal_distribution<> Ndistr(tc_TT, 0.05); // deviation 0.2 * mean
            time = Ndistr(generator);
        }
        else if (LB[tba_i] == 1 && LB[tba_j] == 1)
        {
            std::normal_distribution<> Ndistr(tc_PP, 0.05); // deviation 0.2 * mean
            time = Ndistr(generator);
        }
        else
        {
            std::normal_distribution<> Ndistr(tc_TP, 0.05); // deviation 0.2 * mean
            time = Ndistr(generator);
        }
    }
    Tc = time;
}

// ---------------------------------------------------------------------
void Cadherin::removeC(vector<vector<vector<vector<double>>>> &Sc, vector<vector<int>> &nc, vector<vector<vector<double>>> &Tc, int *idx, int N, int Nc)
{
    int i = idx[0];
    int j = idx[1];
    int k = idx[2];
    vector<vector<int>> nnc(N, vector<int>(N, 0));
    vector<vector<vector<vector<double>>>> nSc(N, vector<vector<vector<double>>>(N, vector<vector<double>>(Nc, vector<double>(2, 0))));
    vector<vector<vector<double>>> nTc(N, vector<vector<double>>(N, vector<double>(Nc, 0)));
    nSc = Sc;
    nnc = nc;
    nTc = Tc;

    // update nnc
    nnc[i][j] = nc[i][j] - 1;
    nnc[j][i] = nc[j][i] - 1;

    // update nSc
    for (int x = k; x < Nc - 1; x++)
    {
        for (int l = 0; l < 2; l++)
        {
            nSc[i][j][x][l] = Sc[i][j][x + 1][l];
            nSc[j][i][x][l] = Sc[j][i][x + 1][l];
        }
    }
    nSc[i][j][Nc - 1][0] = 0;
    nSc[i][j][Nc - 1][1] = 0;
    nSc[j][i][Nc - 1][0] = 0;
    nSc[j][i][Nc - 1][1] = 0;

    // update nTc
    for (int x = k; x < Nc - 1; x++)
    {
        nTc[i][j][x] = Tc[i][j][x + 1];
        nTc[j][i][x] = Tc[j][i][x + 1];
    }
    nTc[i][j][Nc - 1] = 0;
    nTc[j][i][Nc - 1] = 0;

    // substitute back
    Sc = nSc;
    nc = nnc;
    Tc = nTc;
}

// ----------------------------------------------------------------------
void Cadherin::addC(vector<vector<vector<vector<double>>>> &Sc, vector<vector<int>> &nc, vector<vector<vector<double>>> &Tc, int *cidx, double *ccor, double ctime, int N, int Nc)
{
    int i = cidx[0];
    int j = cidx[1];

    vector<vector<int>> nnc(N, vector<int>(N, 0));
    vector<vector<vector<vector<double>>>> nSc(N, vector<vector<vector<double>>>(N, vector<vector<double>>(Nc, vector<double>(2, 0))));
    vector<vector<vector<double>>> nTc(N, vector<vector<double>>(N, vector<double>(Nc, 0)));
    nSc = Sc;
    nnc = nc;
    nTc = Tc;

    // only when there's space for the new c-site
    if (nc[i][j] != Nc)
    {
        int position = nc[i][j];
        nSc[i][j][position][0] = ccor[0];
        nSc[j][i][position][0] = ccor[0];
        nSc[i][j][position][1] = ccor[1];
        nSc[j][i][position][1] = ccor[1];
        nTc[i][j][position] = ctime;
        nTc[j][i][position] = ctime;
        nnc[i][j] = nc[i][j] + 1;
        nnc[j][i] = nc[j][i] + 1;
    }

    // substitute back
    Sc = nSc;
    nc = nnc;
    Tc = nTc;
}

/************************************************************************************
                                  Integrin
************************************************************************************/

void Integrin::gen_i(vector<vector<double>> &X, GIpara gi_para, vector<vector<vector<double>>> &Si, vector<int> &ni, vector<vector<double>> &Ti, vector<int> &LB)
{
    // randomly generate new c-sites from cells X
    std::random_device rd;
    std::mt19937 generator(rd());

    int N = gi_para.N;
    int Ni = gi_para.Ni;
    double r = gi_para.r;
    double Li = gi_para.Li;
    vector<double> argC(N, 0);
    argC = gi_para.argC;
    double ti_PST = gi_para.ti_PST;
    double ti_PSP = gi_para.ti_PSP;

    for (int i = 0; i < N; i++)
    {
        if (X[i][1] >= r && X[i][1] <= r + Li)
        {
            if (argC[i] == -pi / 2 || argC[i] == pi / 2)
            {
                ni[i] = 1;
                Si[i][0][0] = X[i][0];
            }
            else
            {
                //   int NNi = 0;
                //   if (LB[i] == 0){
                //       NNi = Ni;
                //   }
                //   else if (LB[i] == 1){
                //       NNi = Ni / 2; // PSP cells have less i-sites?
                //   }
                for (int k = 0; k < Ni; k++)
                {
                    std::uniform_real_distribution<> Udistr(0, pi / 2); // pi/3 --> pi/2 on Jul11
                    double theta = Udistr(generator);                   // angle ~ U(0,pi/4)
                    double len_i = X[i][1] / cos(theta);
                    if (len_i <= r + Li)
                    {
                        ni[i] = ni[i] + 1; // must form an i-site
                        int pos = ni[i] - 1;
                        if (argC[i] > -pi / 2 && argC[i] < pi / 2)
                        { // point to the right
                            Si[i][pos][0] = X[i][0] + tan(theta) * X[i][1];
                        }
                        else
                        {
                            Si[i][pos][0] = X[i][0] - tan(theta) * X[i][1];
                        }
                    }
                }
            }
        }
    }

    // attach time
    std::normal_distribution<> Ndistr_T(ti_PST, 0.05);
    std::normal_distribution<> Ndistr_P(ti_PSP, 0.05);

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < ni[i]; k++)
        {
            double time = 0;
            while (time <= 0)
            {
                if (LB[i] == 0)
                {
                    time = Ndistr_T(generator);
                }
                else if (LB[i] == 1)
                {
                    time = Ndistr_P(generator);
                }
            }
            Ti[i][k] = time;
            if (Ti[i][k] == 0)
            {
                cout << "error: Ti is zero"
                     << "\n";
            }
        }
    }
}

// ----------------------------------------------------------------------
void Integrin::gen_i_init(vector<vector<double>> &X, GIpara gi_para, vector<vector<vector<double>>> &Si, vector<int> &ni, vector<vector<double>> &Ti, vector<int> &LB)
{
    // randomly generate new c-sites from cells X
    std::random_device rd;
    std::mt19937 generator(rd());

    int N = gi_para.N;
    int Ni = gi_para.Ni;
    Ni = Ni / 2;
    double r = gi_para.r;
    double Li = gi_para.Li;
    vector<double> argC(N, 0);
    argC = gi_para.argC;
    double ti_PST = gi_para.ti_PST;
    double ti_PSP = gi_para.ti_PSP;

    for (int i = 0; i < N; i++)
    {
        if (X[i][1] >= r && X[i][1] <= r + Li)
        {
            if (argC[i] == -pi / 2 || argC[i] == pi / 2)
            {
                ni[i] = 1;
                Si[i][0][0] = X[i][0];
            }
            else
            {
                //   int NNi = 0;
                //   if (LB[i] == 0){
                //       NNi = Ni;
                //   }
                //   else if (LB[i] == 1){
                //       NNi = Ni / 2;
                //   }
                for (int k = 0; k < Ni; k++)
                {
                    std::uniform_real_distribution<> Udistr(0, pi / 2); // pi/3 --> pi/2 on Jul11
                    double theta = Udistr(generator);                   // angle ~ U(0,pi/4)
                    double len_i = X[i][1] / cos(theta);
                    if (len_i <= r + Li)
                    {
                        ni[i] = ni[i] + 1; // must form an i-site
                        int pos = ni[i] - 1;
                        if (argC[i] > -pi / 2 && argC[i] < pi / 2)
                        { // point to the right
                            Si[i][pos][0] = X[i][0] + tan(theta) * X[i][1];
                        }
                        else
                        {
                            Si[i][pos][0] = X[i][0] - tan(theta) * X[i][1];
                        }
                    }
                }
            }
        }
    }

    // attach time
    std::uniform_real_distribution<> Udistr_T(ti_PST / 2, ti_PST * 2);
    std::uniform_real_distribution<> Udistr_P(ti_PSP / 2, ti_PSP * 2);

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < ni[i]; k++)
        {
            double time = 0;
            while (time <= 0)
            {
                if (LB[i] == 0)
                {
                    time = Udistr_T(generator);
                }
                else if (LB[i] == 1)
                {
                    time = Udistr_P(generator);
                }
            }
            Ti[i][k] = time;
            if (Ti[i][k] == 0)
            {
                cout << "error: Ti is zero"
                     << "\n";
            }
        }
    }
}

// ----------------------------------------------------------------------
void Integrin::gen_i_single(vector<vector<double>> &X, int tba_i, GIpara gi_para, double *Si, double &Ti, vector<int> &LB)
{
    // intg.gen_i_single(X, tba_i, gi_para, icor_add, itime_add, LB);
    std::random_device rd;
    std::mt19937 generator(rd());

    int N = gi_para.N;
    int Ni = gi_para.Ni;
    double r = gi_para.r;
    double Li = gi_para.Li;
    vector<double> argC(N, 0);
    argC = gi_para.argC;
    double ti_PST = gi_para.ti_PST;
    double ti_PSP = gi_para.ti_PSP;
    int i = tba_i;

    if (X[i][1] >= r && X[i][1] <= r + Li)
    {
        if (argC[i] == -pi / 2 || argC[i] == pi / 2)
        {
            Si[0] = X[i][0];
        }
        else
        {
            std::uniform_real_distribution<> Udistr(0, pi / 2); // pi/3 --> pi/2 on Jul11
            double theta = Udistr(generator);
            double len_i = X[i][1] / cos(theta);
            if (len_i <= r + Li)
            {
                if (argC[i] > -pi / 2 && argC[i] < pi / 2)
                { // point to the right
                    Si[0] = X[i][0] + tan(theta) * X[i][1];
                }
                else
                {
                    Si[0] = X[i][0] - tan(theta) * X[i][1];
                }
            }
        }
    }

    // attach time
    double time = 0;
    std::normal_distribution<> Ndistr_T(ti_PST, 0.05); // deviation 0.2 * mean
    std::normal_distribution<> Ndistr_P(ti_PSP, 0.05); // deviation 0.2 * mean

    while (time <= 0)
    {
        if (LB[i] == 0)
        {
            time = Ndistr_T(generator);
        }
        else if (LB[i] == 1)
        {
            time = Ndistr_P(generator);
        }
    }
    Ti = time;
    if (Ti == 0)
    {
        cout << "error: Ti is zero"
             << "\n";
    }
}

// ---------------------------------------------------------------------
void Integrin::removeI(vector<vector<vector<double>>> &Si, vector<int> &ni, vector<vector<double>> &Ti, int *idx, int N, int Ni)
{
    int i = idx[0];
    int k = idx[1];
    vector<int> nni(N, 0);
    vector<vector<vector<double>>> nSi(N, vector<vector<double>>(Ni, vector<double>(2, 0)));
    vector<vector<double>> nTi(N, vector<double>(Ni, 0));
    nSi = Si;
    nni = ni;
    nTi = Ti;

    // update nni
    nni[i] = ni[i] - 1;

    // update nSi
    for (int x = k; x < Ni - 1; x++)
    {
        for (int l = 0; l < 2; l++)
        {
            nSi[i][x][l] = Si[i][x + 1][l];
        }
    }
    nSi[i][Ni - 1][0] = 0;
    nSi[i][Ni - 1][1] = 0;

    // update nTi
    for (int x = k; x < Ni - 1; x++)
    {
        nTi[i][x] = Ti[i][x + 1];
    }
    nTi[i][Ni - 1] = 0;

    // substitute back
    Si = nSi;
    ni = nni;
    Ti = nTi;
}

// ----------------------------------------------------------------------
void Integrin::addI(vector<vector<vector<double>>> &Si, vector<int> &ni, vector<vector<double>> &Ti, int idx, double *icor, double itime, int N, int Ni)
{
    int i = idx;

    vector<int> nni(N, 0);
    vector<vector<vector<double>>> nSi(N, vector<vector<double>>(Ni, vector<double>(2, 0)));
    vector<vector<double>> nTi(N, vector<double>(Ni, 0));
    nSi = Si;
    nni = ni;
    nTi = Ti;

    // only when there's space for the new c-site
    if (ni[i] != Ni)
    {
        int pos = ni[i];
        nSi[i][pos][0] = icor[0];
        nSi[i][pos][1] = icor[1];
        nTi[i][pos] = itime;
        nni[i] = ni[i] + 1;
    }

    // substitute back
    Si = nSi;
    ni = nni;
    Ti = nTi;
}