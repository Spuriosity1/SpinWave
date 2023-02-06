/**
 * @file json_parse_test.cpp
 * @author Alaric Sanders (you@domain.com)
 * @brief Tests correct loading and storage of json lattice specifications
 * @version 0.1
 * @date 2023-02-06
 * 
 * @copyright GPLv3
 * 
 */

#include <lattice.hpp>
#include <iostream>

int main(int argc, char** argv){
    lattice L;
    lattice L2;
    if (argc < 3) {
        throw std::runtime_error("Usage: json_parse_test [infile] [outfile]");
    }
    std::string source(argv[1]);
    std::string dest(argv[2]);
    load_json(source, L );

    std::cout.width(4);
    std::cout << "###################################################\n";
    std::cout << "Loaded following from file:\n";
    std::cout << L << "\n\n";
    
    save_json(dest,   L );
    load_json(dest,   L2);
    std::cout<<"###################################################";
    std::cout << "\nReloaded following from file:\n";
    std::cout<< L << "\n\n";

    int res = ( L == L2 );
    std::cout << "Saved object "<<(res==0 ? "is identical to" : "differs from") << " object loaded from file" << std::endl;
    return res;
}
