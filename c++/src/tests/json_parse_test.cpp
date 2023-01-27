#include <lattice.h>

int main(int argc, char** argv){
    lat::lattice L;
    std::string s(argv[1]);
    std::string s2(argv[2]);
    L.load_json(s);
    L.save_json(s2);
}
