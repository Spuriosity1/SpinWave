#ifndef bond_hh
#define bond_hh


#include <coupling_type.hpp>
#include <armadillo>


#ifndef MACHINE_EPS
#define MACHINE_EPS 1e-16
#endif

struct bond;
struct site;

// Represents the name and position of a site
struct site {
    site(size_t i, const std::string& name):idx(i), name(name){}
    // site index
    const size_t idx;
    // A (unique) name for the site
    const std::string name;
    // The position of the site *in lattice coordinates*
    arma::dvec3 xyz;
    // Heisenberg spin of the site
    arma::dvec3 spin;
    // local background magnetic field
    arma::dvec3 bg_H;
    // Neighbours
    // bonds in which this object is 'to'
    std::vector<std::shared_ptr<bond> > neighbours_as_to;
    // bonds in which this object is 'from'
    std::vector<std::shared_ptr<bond> > neighbours_as_from;
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
    coupling(b.coupling), from_site(b.from_site), to_site(b.to_site),
    dx(b.dx)
    {}
    bond(const bond&& b) : 
    coupling(b.coupling), from_site(b.from_site), to_site(b.to_site),
    dx(std::move(b.dx))
    {}
    /**
     * @brief Construct a new bond object.
     * 'from' and 'to' are defined such that 
     * dx = R_to - R_from + [unit cells]
     * 
     * @param coupling Const reference to a coupling, which will be watched.
     * @param from_site Site that the bond is understood to come 'from'
     * @param to_site Site that the bond is understood to go 'to'
     * @param to_cell Unit cell that the bond ends in, if it is thought to start in the cell (0,0,0).
     */
    bond(const coupling_type& coupling, const site& from_site, const site& to_site, const arma::ivec3& to_cell) :
        coupling(coupling), from_site(from_site), to_site(to_site), dx(to_site.xyz + to_cell - from_site.xyz)
        {};

    bond& operator=(const bond&) = delete; // copy assignment forbidden

    const coupling_type& coupling;
    // Sparse array of coupling vectors
    // Interpret these as (site_idx) (vectors R)
    // const size_t from_idx;
    // const size_t to_idx;

    const site& from_site;
    const site& to_site;
    // vector [R_from - R_to]
    arma::dvec3 dx;

    // physics: FROM . M . TO
    double classical_energy() const {
        return arma::dot(from_site.spin, coupling.mat * to_site.spin);
    }
};

inline bool operator==(const bond& b1, const bond& b2) = delete;

inline bool operator!=(const bond& b1, const bond& b2) = delete;

inline std::ostream& operator<<(std::ostream& o, const bond& b){
    const coupling_type& c = b.coupling;
    o << c.name;
    o << " " << b.from_site.name << "->" << b.from_site.name;
    o <<", dx="<<b.dx.t();
    return o;
}

#endif