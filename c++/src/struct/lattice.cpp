/* Responsible for storing the gemetry of a unit cell.
*/

#include "lattice.h"
#include <fstream>

using namespace lat;

template<typename T>
inline arma::vec3 as_armavec(const T& x){
    return arma::vec3({x[0],x[1],x[2]});
}

// RECIPROCALS

const double PI = 3.141592653589793238462643383279502884L;

// returns b(i, n) such that a(m,i)b(i,n) = 2pi \delta_nm
void lattice::reciprocal_vectors(arma::dmat33& b, const arma::dmat33& a){
    double vol = arma::det(a);
    if (fabs(vol) < 1e-16){
        throw "Coordiante system is singular";
    }
    b.col(0) = 2*PI*arma::cross(a.row(2),a.row(3))/vol;
    b.col(1) = 2*PI*arma::cross(a.row(3),a.row(1))/vol;
    b.col(2) = 2*PI*arma::cross(a.row(1),a.row(2))/vol;
}


// CONSTRUCTION
uint16_t lattice::add_site(const std::string& name, const arma::vec3& where){
    site s;
    s.name = name;
    s.xyz = where;
    size_t idx = sites.size();
    this->sites.push_back(s);
    auto res = this->siteDict.insert(std::pair<std::string, size_t>(s.name, idx));
    if (res.second == false) {
        // have a duplicate!
        throw "Attempting to add a duplicate site to map!";
    }
}

// sets the coupling matrix of "handle" to "J"
// passed by reference, so couplings can still be accessed
void lattice::set_coupling(const coupling_type& J){
    this->coupling_types.push_back(J);
    this->couplingDict[J.name] = this->coupling_types.back();
}



void add_bond(size_t from_idx, size_t to_idx, const arma::vec3& to_cell, const coupling_type& J){
    bond b(&J);
    b.from_idx = from_idx, b.to_idx = to_idx;
    throw std::out_of_range("Could not add bond to site");
}

void add_bond(const arma::dvec3& from_loc, const arma::dvec3& to_loc){
    throw std::out_of_range("Could not add bond to site");
}



void lattice::clear(){
    sites.clear();
    siteDict.clear();
    coupling_types.clear();
    bonds.clear();
}


///// INPUT/OUTPUT

bool lattice::load_json(const std::string& filename){
    std::ifstream ifs(filename);
    json jf = json::parse(ifs, nullptr, true, true); // parse ignoring comments

    // delete everything
    clear();

    // Cursed parser code
    try {
        // Pull out the subarrays
        copy_to_mat(this->lattice_vectors, jf.at("lattice vectors"));
        copy_to_mat(this->magnetic_vectors, jf.at("magnetic vectors"));
            
        // populate sites
        for (auto& s : jf.at("sites")) {
            add_site(s.at("name"), as_armavec<json::array_t>(s.at("xyz")));
        }

        // define couplings
        for (auto&& [name, ctype]: jf.at("couplings").items()) {
            coupling_type J;
            J.name = name;
            J.val = 0;

            json::array_t arrarr = ctype.at("matrix");
            copy_to_mat(J.mat, arrarr);

            this->coupling_types.push_back(J);
            // populate bond list
            for (auto&& bondspec : ctype.at("bonds") ){
                bond b(&J);
                b.from_idx = siteDict[bondspec["from"]];
                b.to_idx = siteDict[bondspec["to"]];
                b.dx = sites[b.to_idx].xyz - sites[b.from_idx].xyz;
                arma::dvec3 delta = {bondspec["celldelta"][0], bondspec["celldelta"][1],bondspec["celldelta"][2]};
                b.dx += delta*lattice_vectors;
                this->bonds.push_back(b);
            }    
        }

    } catch (json::out_of_range& e) {
         std::cerr << "JSON file" << filename << "has invalid structure:\n" << e.what()<<std::endl;
         return false;
    }
    return true;
}

void lattice::save_json(const std::string& filename) {
    std::ofstream ofs;
    json j;
    j["lattice vectors"] = json::array();
    j["magnetic vectors"] = json::array();
    for (size_t i=0; i<3; i++){
        j["lattice vectors"][i]  = json::array({lattice_vectors(i,0), lattice_vectors(i,1), lattice_vectors(i,2)});
        j["magnetic vectors"][i] = json::array({magnetic_vectors(i,0),magnetic_vectors(i,1),magnetic_vectors(i,2)});
    }
    j["sites"] = json::array();
    for (size_t i=0; i<this->sites.size(); i++){
        j["atoms"][i] = json({});
        j["atoms"][i]["name"] = sites[i].name;
        j["atoms"][i]["xyz"] = json::array({sites[i].xyz[0],sites[i].xyz[1],sites[i].xyz[2]});   
        j["atoms"][i]["gs_vector"] = sites[i].spin;
    }

    j["couplings"] = json::object();
    for (auto&c : this->coupling_types){
        j["couplings"][c.name] = json::object();
        j["couplings"][c.name]["matrix"] = c.mat;
        j["couplings"][c.name]["bonds"] = json::array();
    }

    for (auto&b : this->bonds){
        auto tmp = json::object();
        
        tmp["from"] = sites[b.from_idx].name;
        tmp["to"] = sites[b.to_idx].name;
        // dx = to - from + celldelta*latvecs
        // celldelta*latvecs = (dx -to + from)
        // latvecsT *celldeltaT = (dc-to+from)T
        // determine the unit cell
        arma::dvec3 v = arma::solve(lattice_vectors.t(), sites[b.from_idx].xyz - sites[b.to_idx].xyz + b.dx);
        tmp["celldelta"] =  json::array({lround(v[0]), lround(v[1]), lround(v[2])}) ;

       j["couplings"][b.coupling->name]["bonds"].push_back(tmp);
    }

    ofs << j;
    ofs.close();
}

// utility
// copies json array 
void lattice::copy_to_mat(arma::dmat33& m, const json::array_t& a){
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            m(i,j) = a[i][j];
        }
    }
}


