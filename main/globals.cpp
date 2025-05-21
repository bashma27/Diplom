#include "globals.h"

vector<vector<Coord2>> zone_perf;
vector<pair<int, vector<double>>> W;

int NUM_ZONE_PERF = 0;
int NUM_SPLIT_X = 0, NUM_SPLIT_Y = 0, NUM_SPLIT_Z = 0;
int NUM_NODES_IN_EDGE_X = 0, NUM_NODES_IN_EDGE_Y = 0, NUM_NODES_IN_EDGE_Z = 0, NUM_NODES = 0;
int num_ph = 0;

double EPS = 1e-300, F = 0;

vector<vector<double>> G = {
    {7., -8., 1.},
    {-8., 16., -8.},
    {1., -8., 7.}
};
vector<vector<double>> M = {
    {4., 2., -1.},
    {2., 16., 2.},
    {-1., 2., 4.}
};

vector<double> set_flow_zp;
vector<int> ia, ja, choice;
vector<double> aal, di, b, q, L_sq, di_sq, normal;
vector<vector<double>> k_ph;
vector<double> K, eta_ph;
vector<Coord3> nodes;
unordered_set<int> face_1;
vector<pair<int, vector<int>>> face_2_zp;
vector<pair<int, vector<int>>> array_p;
vector<int> ident_fict;
vector<int> i_ident_fict;
vector<int> fict_nodes;
vector<double> faces_flow_value;
vector<vector<double>> faces_flow_ph_value;
vector<pair<int, int>> num_faces_zp;
vector<int> fict_faces;
vector<int> fict_el;
vector<double> q_flow;
vector<vector<int>> list_face;
vector<int> el_for_face;
vector<vector<double>> satur;