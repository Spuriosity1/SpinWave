/**
 * @file luttinger.hpp
 * @author your name (you@domain.com)
 * @brief Luttinger-Tisza optimisation of the classical ground state.
 * @version 0.1
 * @date 2023-02-06
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef OPTIMISE_CXX_H
#define OPTIMISE_CXX_H


#include <lattice.hpp>
#include <string>
#include <random>

struct optimise_result {
    std::string msg;
    double energy_per_cell;
    double var_E;
};

double classical_energy(const lattice& L){
    double E;
    for (const bond& b : L.get_bonds()){
        E += b.classical_energy();
    }
    return E;
}

inline arma::dvec3 loc_field (const site& s){
    // H = from * M * to
    arma::dvec3 h_loc;
    for (auto b : s.neighbours_as_from){
        h_loc += b->coupling.mat * b->to_site.xyz;
    }
    for (auto b : s.neighbours_as_to){
        h_loc += b->coupling.mat.t() * b->from_site.xyz;
    }
    return h_loc + s.bg_H;
}

/**
 * @brief Luttinger-Tisza optimisation of classical magnetic moments.
 * 
 * @param L lattice. 'Heisenberg' (i.e. classical) moments are reoriented to a classical ground state.
 * @param optimise_result a struct containing information about the success (or failue) of the optimisation, as well as estimates on the energy error.
 */
void LT_optimise(lattice& L, optimise_result& res){

}

/**
 * @brief Monte Carlo optimisation of classical magnetic moments.
 * 
 * @param L lattice. 'Heisenberg' (i.e. classical) moments are reoriented to a classical ground state.
 * @param optimise_result a struct containing information about the success (or failue) of the optimisation, as well as estimates on the energy error.
 */
void CMC_step(lattice& L, optimise_result& res, std::mt19937& random_engine, double kT){
    // Strategy: small fluctuations guided by gradient descent
    // double E = classical_energy(L);

    std::uniform_real_distribution<double> unif(0, 1);

    arma::dvec3 tmp_h, tmp_s;

    if (kT > 1e-10){
        for (size_t j =0; j<L.num_sites(); j++){
            size_t idx = random_engine() % L.num_sites();
            // Propose the move:
            tmp_h = loc_field(L.get_site(idx));
            tmp_s = normalise(tmp_h + arma::randn(3, 1, arma::distr_param(0.,kT)) );
            if ( unif(random_engine) < exp(-(arma::dot(tmp_h, tmp_s) - arma::dot(tmp_h, L.get_site(idx).spin))/kT)){
                // accept
                L.cl_spin_at(idx) = tmp_s;
            }
        }

    } else {
        // basically zero temperature, pure gradient descent!
        for (size_t j =0; j<L.num_sites(); j++){
            size_t idx = random_engine() % L.num_sites();
            // Propose the move:
            tmp_h = loc_field(L.get_site(idx));
            tmp_s = normalise(tmp_h + arma::randn(3, 1, arma::distr_param(0.,kT)) );
            L.cl_spin_at(idx) = tmp_s;
        }
    }

    
}


void randomise_sites(lattice& L) {
    // randomise directions
    for (size_t j=0; j< L.num_sites(); j++){
        L.cl_spin_at(j) = normalise( arma::randn(3, 1) );
    }
}



#endif