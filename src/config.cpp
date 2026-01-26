#include "config.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;

// Ce fichier implémente la lecture des paramètres de configuration depuis un fichier texte

void SimulationConfig::readFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Impossible d'ouvrir le fichier de configuration: " + filename);
    }
    string line;
    while (getline(file, line)) {
        // On ignore les commentaires et lignes vides
        line = line.substr(0, line.find('#'));
        if (line.empty()) continue;
        istringstream iss(line);
        string key;
        if (getline(iss, key, '=')) {
            string value;
            if (getline(iss, value)) {
                key.erase(remove_if(key.begin(), key.end(), ::isspace), key.end());
                value.erase(remove_if(value.begin(), value.end(), ::isspace), value.end());
                // On utilise un ensemble de conditions pour assigner les valeurs lues aux attributs correspondants
                if (key == "solver_type") {
                    if (value == "NEUTRONIC") solver_type = SolverType::NEUTRONIC;
                    else if (value == "COUPLED") solver_type = SolverType::COUPLED;
                } 
                else if (key == "guide_positions") {
                    guide_positions.clear();
                    istringstream pos_stream(value);
                    string pos_pair;
                    while (getline(pos_stream, pos_pair, ';')) {
                        size_t comma_pos = pos_pair.find(',');
                        if (comma_pos != string::npos) {
                            int x = stoi(pos_pair.substr(0, comma_pos));
                            int y = stoi(pos_pair.substr(comma_pos + 1));
                            guide_positions.push_back(make_pair(x, y));
                        }
                    }
                }
                else if (key == "scheme_type") {
                    if (value == "IMPLICIT_EULER") scheme_type = SchemeType::IMPLICIT_EULER;
                    else if (value == "CRANK_NICOLSON") scheme_type = SchemeType::CRANK_NICOLSON;
                }

                else if (key == "boundary_condition") {
                    if (value == "DIRICHLET") bc_type = BCType::DIRICHLET;
                    else if (value == "NEUMANN") bc_type = BCType::NEUMANN;
                    else if (value == "ROBIN") bc_type = BCType::ROBIN;
                }

                else if (key == "dimension") dimension = stoi(value);
                else if (key == "N") N = stoi(value);
                else if (key == "pitch") pitch = stod(value);
                else if (key == "r_fuel") r_fuel = stod(value);
                else if (key == "r_guide") r_guide = stod(value);
                else if (key == "cells_per_pitch") cells_per_pitch = stoi(value);
                else if (key == "f_ins") f_ins = stod(value);
                else if (key == "tf") tf = stod(value);       
                else if (key == "dt") dt = stod(value);
                else if (key == "n_out") n_out = stoi(value);
                
                else if (key == "initial_phi") initial_phi = stod(value);
                else if (key == "initial_T") initial_T = stod(value);
                else if (key == "length_unit") length_unit = value;
            }
        }
    }
    file.close();
} 

void SimulationConfig::printSummary(const std::string& filename) const {
    using std::cout; using std::endl;
    cout << "Paramètres de la simulation :" << endl;
    cout << " Fichier de configuration : " << filename << endl;
    cout << " Dimensions de l'assemblage : " << N << "x" << N << endl;
    cout << " Nombre de cellules par pitch : " << cells_per_pitch << endl;
    cout << " Facteur d'insertion des barres de contrôle : " << f_ins << endl;
    cout << " Durée totale de la simulation : " << tf << " s" << endl;
    cout << " Pas de temps : " << dt << " s" << endl;
    cout << " Schéma temporel : Implicit Euler" << endl;
}
