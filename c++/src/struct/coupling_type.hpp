
#ifndef coupling_type_hh
#define coupling_type_hh

#include <iostream>
#include <armadillo>

// Defines a type of fundamental interaction (e.g. Heisenberg)
struct coupling_type {
    // coupling_type(const std::string& name, const arma::dmat33& mat) : name(name), mat(mat){};
    std::string name;
    arma::dmat33 mat;
    double val;

};

inline std::ostream& operator<<(std::ostream& o, const coupling_type& c){
    o << "["<<c.name<<"] = "<<c.val<<" *\n"<<c.mat;
    return o;
}

namespace coupling_matrices {

static const arma::dmat33 heis                 = arma::dmat33({{1,0,0},{0,1,0},{0,0,1}});
static const arma::dmat33 ising_z              = arma::dmat33({{0,0,0},{0,0,0},{0,0,1}});
static const arma::dmat33 gamma_z              = arma::dmat33({{0,1,0},{-1,0,0},{0,0,0}});
static const arma::dmat33 gammap_z             = arma::dmat33({{0,0,1},{0,0,1},{1,1,0}});
static const arma::dmat33 xxyy                 = arma::dmat33({{1,0,0},{0,1,0},{0,0,0}});

};



#endif