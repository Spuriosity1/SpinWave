/**
 * @file lattice.hpp
 * @author Alaric Sanders (you@domain.com)
 * @brief Stores the abstractions responsible for representing a lattice's geometry.
 * @version 0.2
 * @date 2023-02-06
 * 
 * @copyright GPLv3
 * 
 */


#ifndef lattice_hh
#define lattice_hh

#include <armadillo>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <nlohmann/json.hpp>

#include <coupling_type.hpp>
#include <bond.hpp>



// namespace lat {

using json=nlohmann::json;




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

    inline size_t num_sites() const {
        return sites.size();
    }

    inline const site& get_site(size_t i) const{
        return sites[i];
    }

    inline const site& get_site(const std::string& s) const{
        return sites[site_dict.at(s)];
    }

    inline arma::dvec3& cl_spin_at(size_t i) {
        return sites[i].spin;
    }

    inline arma::dvec3& cl_spin_at(const std::string& s){
        return sites[site_dict.at(s)].spin;
    }

    bool has_site(const std::string& name) const {
        for (auto& s : sites){
            if (s.name == name) return true;
        }
        return false;
    }

    bool has_coupling(const std::string& name) const {
        for (auto& s : coupling_types){
            if (s.name == name) return true;
        }
        return false;
    }

    inline const coupling_type& get_coupling(const std::string& s) const{
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

    
    

protected:
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
