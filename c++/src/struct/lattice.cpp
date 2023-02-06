/**
 * @file lattice.cpp
 * @author Alaric Sanders (you@domain.com)
 * @brief Stores the abstractions responsible for representing a lattice's geometry.
 * @version 0.2
 * @date 2023-02-06
 * 
 * @copyright GPLv3
 * 
 */

#include "lattice.hpp"
#include <fstream>


using json=nlohmann::json;
// utility

template<typename T, typename vt=arma::dvec3>
inline vt as_armavec(const T& x){
    return vt({x[0],x[1],x[2]});
}




// copies json array 
inline arma::dmat33 dmat33_from_jsonarr(const json::array_t& a){
    arma::dmat33 m;
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            m(i,j) = a[i][j];
        }
    }
    return m;
}

// copies json array 
inline json::array_t jsonarr_from_dmat33(const arma::dmat33& m){
    json::array_t a;
    a.reserve(3);
    for (int i=0; i<3; i++){
        a.push_back(json::array({m(i,0), m(i,1), m(i,2)}));
    }
    return a;
}



/**
 * @brief returns b(i, n) such that a(m,i)b(i,n) = 2pi \delta_nm
 * 
 * @param b Matrix to store reciprocal vectors in
 * @param a Knonw unit  cell dimensions
 */
void reciprocal_vectors(arma::dmat33& b, const arma::dmat33& a) {
    const static double PI = 3.141592653589793238462643383279502884L;
    double vol = arma::det(a);
    if (fabs(vol) < 1e-16){
        throw "Coordiante system is singular";
    }
    b.col(0) = 2*PI*arma::cross(a.row(2),a.row(3))/vol;
    b.col(1) = 2*PI*arma::cross(a.row(3),a.row(1))/vol;
    b.col(2) = 2*PI*arma::cross(a.row(1),a.row(2))/vol;
}



bool operator==(const lattice& l1, const lattice& l2){
    try {
        // Check that site dicts agree on keys
        for (auto const& [name, site_idx]: l1.site_dict) {
            if (l2.get_site(name) != l1.get_site(name)) return false;
        }

        for (auto const& [name, J_idx]: l1.coupling_dict) {
            if (l1.get_coupling(name) != l2.get_coupling(name)) return false;
        }

        if (l1.bonds.size() != l2.bonds.size() ) return false;

        std::list<size_t> idx_set(l1.bonds.size());

        for (size_t i=0; i<l1.bonds.size(); i++){
            idx_set.push_back(i);
        }

        for (auto const& b1 : l1.bonds) {
            // may be reshuffled in l2, must do linear search :(
            bool found_one = false;

            // Search the remaining elements
            for (std::list<size_t>::const_iterator it = idx_set.begin(); it != idx_set.end(); it++){
                auto const& b2 = l2.bonds[*it];
                if (b1.coupling.name == b2.coupling.name &&
                l1.sites[b1.to_idx] == l2.sites[b2.to_idx] &&
                l1.sites[b1.from_idx] == l2.sites[b2.from_idx] &&
                arma::norm(b1.dx - b2.dx, 2) < MACHINE_EPS
                    ) {
                    // found one!
                    found_one = true;
                    idx_set.erase(it);
                    it = idx_set.end();
                }
            }
            if (!found_one) return false;
            // ensures only 1/2 N(N+1) equality tests.
        }

    } catch (std::out_of_range& e) {
        std::cerr << e.what();
        return false;
    }
    return true;
}


void lattice::get_reciprocal_lat_vectors(arma::dmat33& b) const {
    reciprocal_vectors(b, this->lattice_vectors);
}

void lattice::get_reciprocal_mag_vectors(arma::dmat33& b) const {
    reciprocal_vectors(b, this->magnetic_vectors);
}

/**
 * @brief Defines a new site for the lattice
 * 
 * @param name Name for the siet (must be unique)
 * @param where Vector3 indicating location of site 
 * @param lattice_coords Boolean flag for whether 'where' is supplied in lattice coordiantes or Cartesian coordinates
 * @return uint16_t index of inserted site
 */
uint16_t lattice::add_site(const std::string& name, const arma::vec3& where, bool lattice_coords=true){
    site s;
    s.name = name;
    s.xyz = lattice_coords ? where : arma::solve(this->lattice_vectors, where);
    
    size_t idx = sites.size();
    this->sites.push_back(s);
    auto res = this->site_dict.insert(std::pair<std::string, size_t>(s.name, idx));
    // res = (iterator pointing to newly inserted memeber, success flag)
    if (res.second == false) {
        // have a duplicate!
        throw std::logic_error("Attempting to add a duplicate site - all sites must be given unique names.");
    }
    return idx;
}


void lattice::add_bond(const std::string &from, const std::string& to, const arma::Col<arma::sword>::fixed<3>& to_cell, const std::string& J_name ){
    const coupling_type& J = get_coupling(J_name);
    try {
        add_bond(this->site_dict.at(from), this->site_dict.at(to), to_cell, J );
    } catch (const std::out_of_range& e) {
        std::ostringstream ss;
        ss <<e.what()<<"\nFailed to add bond: Could not resolve bondspec [" <<from<<"] -> ["<<to<<"], "<<J_name<<std::endl;
        throw std::out_of_range(ss.str());
    }
}

/**
 * @brief Adds a bond linking site [from_idx] in the original unit cell to site [to_idx] in unit cell of index to_cell [hkl] 
 * 
 * @param from_idx Site to link from
 * @param to_idx Site to link to
 * @param to_cell [hkl] index of site to link to
 * @param J Pointer to a coupling_type object
 */
void lattice::add_bond(size_t from_idx, size_t to_idx, const arma::Col<arma::sword>::fixed<3>& to_cell, const coupling_type& J){
    bond b(J);
    b.from_idx = from_idx, b.to_idx = to_idx;
    b.dx = this->sites[from_idx].xyz + to_cell - this->sites[to_idx].xyz;
    // sanity check: Is J one of my couplings?
    for (auto& c: this->coupling_types){
        if (c.name == J.name){
            this->bonds.push_back(b);
            return;
        }
    }
    throw std::logic_error("Coupling type specified does not belong to this lattice!");
}


// void lattice::add_bond(const arma::dvec3& from_loc, const arma::dvec3& to_loc){
//     throw std::out_of_range("Could not add bond to site");
// }

void lattice::define_coupling(coupling_type& J){
    coupling_dict.insert({J.name, coupling_types.size()});
    coupling_types.push_back(J);
}

// void lattice::define_coupling(coupling_type&& J){
//     std::string jn = J.name;
//     coupling_dict.insert(
//         std::pair(jn, coupling_types.size())
//     );
//     coupling_types.push_back(J);
// }

void lattice::clear(){
    sites.clear();
    site_dict.clear();
    bonds.clear();
}

///// INPUT/OUTPUT
bool load_json(const std::string& filename, lattice& lat){
    std::ifstream ifs(filename);
    json jf = json::parse(ifs, nullptr, true, true); // parse ignoring comments

    // delete everything
    lat.clear();

    // Cursed parser code
    // TODO: data validation. At this point it's spray & pray.
    try {
        
        for (int i=0; i<3; i++){
            arma::dvec3 v = as_armavec<json::array_t, arma::dvec3>(jf.at("lattice vectors")[i]);
            lat.set_lattice_vector(i, v);
            arma::dvec3 mv = as_armavec<json::array_t, arma::dvec3>(jf.at("magnetic vectors")[i]);
            lat.set_magnetic_vector(i, mv);
        }

        bool use_lattice_vectors = false;
        if (jf.at("site_basis") == "lattice"){
            use_lattice_vectors = true;
        } else if (jf.at("site_basis") == "cartesian"){
            use_lattice_vectors = false;
        } else {
            throw std::runtime_error("Expected 'lattice' or 'cartesian' for specifier file");
        }

        // populate sites
        for (auto&& s : jf.at("sites")) {
            lat.add_site(s.at("name"), as_armavec<json::array_t>(s.at("xyz")), use_lattice_vectors);
        }

        // define couplings
        for (auto& [name, ctype]: jf.at("couplings").items()) {
            coupling_type J =coupling_type();
            J.name = name;
            J.val = ctype.at("strength");

            const json::array_t& arrarr = ctype.at("matrix");
            J.mat = dmat33_from_jsonarr(arrarr);

            lat.define_coupling(J);
        }

        for (auto& [name, ctype]: jf.at("couplings").items()) {
            // populate bond list    
            for (auto& bondspec : ctype.at("bonds") ){
                lat.add_bond(bondspec["from"], bondspec["to"], 
                    as_armavec<json::array_t, arma::Col<arma::sword>::fixed<3> >(bondspec["celldelta"]), name);
            }
        }

        

    } catch (json::out_of_range& e) {
         std::cerr << "JSON file '" << filename << "' has invalid structure:\n" << e.what() << std::endl;
         return false;
    }
    return true;
}


// serialisers
inline json as_json(const site& s){
    json j;
    j["name"] = s.name;
    j["xyz"] = json::array({s.xyz[0],s.xyz[1],s.xyz[2]});
    j["heis_vector"] = json::array({s.spin[0],s.spin[1],s.spin[2]});
    return j;
}

inline json as_json(const coupling_type& ct){
    json j;
    j["name"] = ct.name;
    j["matrix"] = jsonarr_from_dmat33( ct.mat);
    j["strength"] = ct.val;
    return j;
}

void save_json(const std::string& filename, const lattice& lat) {
    
    json j;
    j["site_basis"] = "lattice";
    j["lattice vectors"] = json::array();
    j["magnetic vectors"] = json::array();
    for (size_t i=0; i<3; i++){
        arma::drowvec3 v;
        v = lat.get_lattice_vector(i);
        j["lattice vectors"][i]  = json::array({v(0), v(1), v(2)});
        v = lat.get_magnetic_vector(i);
        j["magnetic vectors"][i]  = json::array({v(0), v(1), v(2)});
    }
    j["sites"] = json::array();
    for (size_t i=0; i<lat.num_sites(); i++){
        j["sites"][i] = as_json(lat.get_site(i)); 
    }

    j["couplings"] = json::object();
    for (auto& c : lat.get_coupling_types()){
        j["couplings"][c.name] = as_json(c);
        j["couplings"][c.name]["bonds"] = json::array();
    }

    for (auto& b : lat.get_bonds()){
        json::object_t tmp = json::object();

        const site& from = lat.get_site(b.from_idx);
        const site& to = lat.get_site(b.to_idx);
        
        tmp["from"] = from.name;
        tmp["to"] = to.name;
        // dx = to - from + celldelta
        // celldelta*latvecs = (dx -to + from)
        // determine the unit cell
        // internal units are lattice units: no need to do any matrix inversion
        arma::dvec3 v =  from.xyz - to.xyz + b.dx;
        tmp["celldelta"] =  json::array({lround(v[0]), lround(v[1]), lround(v[2])}) ;
        
        j["couplings"][b.coupling.name]["bonds"].push_back(tmp);
    }

    std::ofstream ofs(filename);
    ofs << j;
    ofs.close();
}

std::ostream& operator<<(std::ostream& o, const lattice& lat){
    o<<"Sites:\n";
    for (auto& s: lat.sites){
        o << "\t"<< s << "\n";
    }
    o<<"Bonds:\n";
    for (auto& b: lat.bonds){
        o << "\t" << b <<"\n";
    }
    o<<"Couplings:\n";
    for (auto& c: lat.coupling_types){
        o << "\t" << c <<"\n";
    }
    return o;
}