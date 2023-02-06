#ifndef lattice_hh
#define lattice_hh

#include <armadillo>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <nlohmann/json.hpp>

#ifndef MACHINE_EPS
#define MACHINE_EPS 1e-16
#endif

// namespace lat {

using json=nlohmann::json;

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


// Represents the name and position of a site
struct site {
    // A (unique) name for the site
    std::string name;
    // The position of the site *in lattice coordinates*
    arma::dvec3 xyz;
    // Heisenberg spin of the site
    arma::dvec3 spin;
};

inline bool operator==(const site& s1, const site& s2){
    return (s1.name == s2.name && arma::norm(s1.xyz - s2.xyz, 2) < MACHINE_EPS && arma::norm(s1.spin - s2.spin, 2) < MACHINE_EPS);
}

inline bool operator!=(const site& s1, const site& s2){
    return !(s1 == s2);
}

inline bool operator==(const coupling_type& J1, const coupling_type& J2){
    return (J1.name == J2.name && arma::norm(J1.mat - J2.mat, 2) < MACHINE_EPS && abs(J1.val - J2.val) < MACHINE_EPS);
}

inline bool operator!=(const coupling_type& J1, const coupling_type& J2){
    return !(J1 == J2);
}

inline std::ostream& operator<<(std::ostream& o, const site& s){
    o << "["<<s.name << "] " <<" @ ["<< s.xyz[0] <<"a1 + "<<s.xyz[1]<<"a2 + " << s.xyz[2]<<"a3]  =" << s.spin.t();
    return o;
}



struct bond {
    bond(const bond& b) : 
    coupling(b.coupling), from_idx(b.from_idx), to_idx(b.to_idx), dx(b.dx)
    {}
    bond(const bond&& b) : 
    coupling(b.coupling), from_idx(b.from_idx), to_idx(b.to_idx)
    {dx = std::move(b.dx);}
    bond(const coupling_type& coupling): coupling(coupling){};
    bond(const coupling_type& coupling, size_t from_idx, size_t to_idx, const arma::dvec3& dx) :
        coupling(coupling), from_idx(from_idx), to_idx(to_idx), dx(dx){};

    bond& operator=(const bond&) = delete; // copy assignment forbidden


    const coupling_type& coupling;
    // Sparse array of coupling vectors
    // Interpret these as (site_idx) (vectors R)
    size_t from_idx;
    size_t to_idx;
    // vector [R_from - R_to]
    arma::dvec3 dx;
};

inline bool operator==(const bond& b1, const bond& b2){
    return (b1.coupling.name == b2.coupling.name \
        && b1.to_idx == b2.to_idx && b1.from_idx == b2.from_idx \
        && arma::norm(b1.dx - b2.dx, 2)<MACHINE_EPS);
}

inline bool operator!=(const bond& b1, const bond& b2){
    return !(b1 == b2);
}

inline std::ostream& operator<<(std::ostream& o, const bond& b){
    const coupling_type& c = b.coupling;
    o << c.name;
    o << " " << b.from_idx << "->" << b.to_idx;
    o <<", dx="<<b.dx.t();
    return o;
}



// Responsible for managing cell geometry
class lattice {
public:
    // Constructs a lattice with unit vectors (1,0,0), (0,1,0), (0,0,1)
    lattice() :
        lattice_vectors(arma::eye(3,3)),
        magnetic_vectors(arma::eye(3,3)) {};
    
    // Class constructor 
    // 'a' contians the three primitive lattice vectors, first index chooses the vector, second is xyz 
    lattice(const arma::dmat33& a) :
        lattice_vectors(a),
        magnetic_vectors(a) {};
        
    lattice(const arma::dmat33& a, const arma::dmat33& a_mag) :
        lattice_vectors(a),
        magnetic_vectors(a_mag) {};


    // friends

    friend bool operator==(const lattice& l1, const lattice& l2);
    friend std::ostream& operator<<(std::ostream& o, const lattice& lat);

    
    // // Adds a spin to the unit cell at postion 'r' (returns index)
    uint16_t add_site(const std::string& sitelabel, const arma::vec3& r, bool lattice_coords);

    void add_bond(const std::string &from, const std::string& to, const arma::Col<arma::sword>::fixed<3>& to_cell, const std::string& J_name );
    void add_bond(size_t from_idx, size_t to_idx, const arma::Col<arma::sword>::fixed<3>& to_cell, const coupling_type& J);
    // void add_bond(const arma::dvec3& from_loc, const arma::dvec3& to_loc);

    // // stores reciprocal vectors in b
    // void reciprocal_lat_vectors(arma::dmat33& b) const;
    // void reciprocal_mag_vectors(arma::dmat33& b) const;

    size_t num_sites() const {
        return sites.size();
    }

    const site& get_site(size_t i) const{
        return sites[i];
    }

    const site& get_site(const std::string& s) const{
        return sites[site_dict.at(s)];
    }

    const coupling_type& get_coupling(const std::string& s) const{
        return coupling_types[coupling_dict.at(s)];
    } 

    // delete everything
    void clear();

    // Lattice vector access control
    inline const arma::drowvec3 get_lattice_vector(arma::uword j) const {
        return lattice_vectors.row(j);
    }

    inline const arma::drowvec3 get_magnetic_vector(arma::uword j) const {
        return magnetic_vectors.row(j);
    }

    // TODO return by value
    void get_reciprocal_lat_vectors(arma::dmat33& b) const;
    void get_reciprocal_mag_vectors(arma::dmat33& b) const;

    inline const arma::dmat33 get_lattice_vectors() const {
        return lattice_vectors;
    }

    inline const arma::dmat33 get_magnetic_vectors() const {
        return magnetic_vectors;
    }



    void set_lattice_vector(arma::uword j, const arma::vec3& v) {
        for (arma::uword k=0; k<3; k++){
            lattice_vectors(j, k) = v(k);
        }
    }

    void set_magnetic_vector(arma::uword j, const arma::vec3& v) {
        for (arma::uword k=0; k<3; k++){
            magnetic_vectors(j,k) = v(k);
        }
    }

    

    void define_coupling(coupling_type& J);
    // void define_coupling(coupling_type&& J);

    inline const std::vector<coupling_type>& get_coupling_types() const {
        return this->coupling_types;
    }
    inline const std::vector<bond>& get_bonds() const {
        return this->bonds;
    }
    inline const std::vector<site>& get_sites() const {
        return this->sites;
    }

    
    

private:
    // Lattice vectors, used for defining BZ coordinates
    // the [i]th vector is in lattice_vectors(i,*);
    arma::dmat33 lattice_vectors;
    // Magnetic unit cell vectors, used internally
    arma::dmat33 magnetic_vectors;
    
	// Stores the positions and GS orientations of [first unit cell] sites
	std::vector<site> sites;
    std::vector<coupling_type> coupling_types;
    std::vector<bond> bonds;

    // dicts to help
    std::map<std::string, size_t> site_dict;
    std::map<std::string, size_t> coupling_dict;

};


bool load_json(const std::string& filename, lattice& l);
void save_json(const std::string& filename, const lattice& l);



// // File IO: reads json from file, returns true on success
// bool load_json(const std::string& filename, lattice& l);
// // File IO: writes object to file
// void save_json(const std::string& filename, const lattice& l);

// };// end of namespace

#endif
