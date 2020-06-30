//
//  main.cpp
//  GONEaverage
//
//  Created by Enrique Santiago on 18/11/19.
//  Copyright © 2019 Enrique Santiago. All rights reserved.
//
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#define maxlin 10000

int main(int argc, const char * argv[]) {
    int i, a, nlin=0, repeticiones;
    double b,NeA[maxlin],NeG[maxlin],NeH[maxlin];
    int generacion[maxlin];
    
    if (argc!=4){
        std::cerr << " Usage: "<< argv[0]<< " [directory_of_temp_files] [name_of_input_file] [number_of_repeats] "<<std::endl;
        return 1;}
    
    std::string directorio=argv[1];
    std::string fichero=argv[2];
    repeticiones=std::atoi(argv[3]);
    std::ifstream entrada;
    std::ofstream salida;
    std::string line;
    
    for (i=0;i<maxlin;++i){NeA[i]=0.0;NeG[i]=1.0;NeH[i]=0.0;}
    
    for (i=1;i<=repeticiones;++i){
        std::string s = std::to_string(i);
        std::string fichentra=directorio+"/"+fichero+"_"+s+"_GONE_Nebest";
        entrada.open(fichentra, std::ios::in);
        if (!(entrada.is_open())){
            std::cerr << "  Error opening file "<< fichentra <<"\n";
            return 1;}
        nlin=0;
        while (getline(entrada, line)){
            std::istringstream iss(line);
            if (!(iss >> a >> b)){
                std::cerr << " Format error in file "<< fichentra<<std::endl;
                entrada.close();
                return 1;}
            generacion[nlin]=a;
            NeG[nlin]*=pow(b,double(1.0/repeticiones));
            NeH[nlin]+=1.0/b;
            NeA[nlin]+=b;
            ++nlin;
        }
        entrada.close();
    }
    for (i=0;i<nlin;++i){
        NeH[i]=repeticiones/NeH[i];
        NeA[i]=NeA[i]/repeticiones;
    }
    std::string fichsal="Average_"+fichero;
    salida.open(fichsal, std::ios::out);
    salida << "Ne averages over "<<repeticiones<<" independent estimates."<<"\n";
    salida << "Generation"<<"\t"<<"Arithmetic_mean"<<"\t"<<"Geometric_mean"<<"\t"<<"Harmonic_mean"<<"\n";
    for (i=0;i<nlin;++i){
        salida << generacion[i]<<"\t"<<NeA[i]<<"\t"<<NeG[i]<<"\t"<<NeH[i]<<"\n";
    }
    salida.close();
    
    
    return 0;
}

