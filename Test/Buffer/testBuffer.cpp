
#include<iostream>
#include "../../Buffer.hpp"

int main(){
    std::string name("SPHERE");
    std::array<double,4> array4 {1,2,3,4};
    Buffer mybuffer;
    mybuffer.pack(name);
    mybuffer.pack(array4);
    name.clear();
    array4.fill(0);

    mybuffer.unpack(name);
    mybuffer.unpack(array4);
    std::cout<<name<< " "<<array4[0]<<array4[1]<<array4[2]<<array4[3]<<std::endl;
    return 0;    
}