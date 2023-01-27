#ifndef lattice_hh
#define lattice_hh

#include <armadillo>
#include <string>
#include <vector>
#include <map>
#include <nlohmann/json.hpp>


namespace lat {

using json=nlohmann::json;

// Defines a type of fundamental interaction (e.g. Heisenberg)
struct coupling_type {
    // coupling_type(const std::string& name, const arma::dmat33& mat) : name(name), mat(mat){};
    std::string name;
    arma::dmat33 mat;
    double val;
};


// Represents the name and position of a site
struct site {
    std::string name;
    arma::dvec3 xyz;
    arma::dvec3 spin;
};


struct bond {
    bond(const coupling_type* coupling):coupling(coupling){};
    bond(const coupling_type* coupling, size_t from_idx, size_t to_idx, const arma::vec3& dx) :
        coupling(coupling), from_idx(from_idx), to_idx(to_idx), dx(dx){};


    const coupling_type* coupling;
    // Sparse array of coupling vectors
    // Interpret these as (site_idx) (vectors R)
    size_t from_idx;
    size_t to_idx;
    arma::dvec3 dx;
};

// serialisers
json as_json(site s);
json as_json(coupling_type ct);
json as_json(bond b);

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
    
    // // Adds a spin to the unit cell at postion 'r' (returns index)
    uint16_t add_site(const std::string& sitelabel, const arma::vec3& r );

    // // sets the coupling matrix of "handle" to "J"
    void set_coupling(const coupling_type& J);

    void add_bond(size_t from_idx, size_t to_idx, const arma::vec3& to_cell, const coupling_type& J);
    void add_bond(const arma::dvec3& from_loc, const arma::dvec3& to_loc);

    // stores reciprocal vectors in b
    void reciprocal_lat_vectors(arma::dmat33& b) {
        reciprocal_vectors(b, this->lattice_vectors);
    }
    void reciprocal_mag_vectors(arma::dmat33& b) {
        reciprocal_vectors(b, this->magnetic_vectors);
    }

    // File IO: reads json from file, returns true on success
    bool load_json(const std::string& filename);
    // File IO: writes object to file
    void save_json(const std::string& filename);

    size_t num_sites(){
        return sites.size();
    }

    // delete everything
    void clear();

private:
    // Lattice vectors, used for defining BZ coordinates
    // the [i]th vector is in lattice_vectors(i,*);
    arma::dmat33 lattice_vectors;
    // Magnetic unit cell vectors, used internally
    arma::dmat33 magnetic_vectors;
    
	// Stores the positions and GS orientations of [first unit cell] sites
	std::vector<site> sites;
    std::map<std::string, size_t> siteDict;


    std::vector<coupling_type> coupling_types;
    std::map<std::string, coupling_type&> couplingDict;

    std::vector<bond> bonds;

    // utility function
    static void copy_to_mat(arma::dmat33& m, const json::array_t& a);

    void reciprocal_vectors(arma::dmat33& b, const arma::dmat33& a);


    // returns the unit cell that x lives in.
    arma::vec3 cell_of(arma::dvec3 x);

};

};// end of namespace

#endif
