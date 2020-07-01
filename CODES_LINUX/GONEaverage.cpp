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
    int i, a, nlin=0, nlin2=0, repeticiones;
    double b,c,d2o,d2p,NeA[maxlin],NeG[maxlin],NeH[maxlin];
    double dif2=0.0,crec[maxlin],d2obs[maxlin],d2prd[maxlin];
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
    
    for (i=0;i<maxlin;++i){
        NeA[i]=0.0;NeG[i]=1.0;NeH[i]=0.0;
        crec[i]=0;d2obs[i]=0;d2prd[i]=1.0;}
    
    for (i=1;i<=repeticiones;++i){
        std::string s = std::to_string(i);
        // Ne:
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
            if (nlin>=maxlin){break;}
        }
        entrada.close();
        // d2:
        std::string fichentra2=directorio+"/"+fichero+"_"+s+"_GONE_d2";
        entrada.open(fichentra2, std::ios::in);
        if (!(entrada.is_open())){
            std::cerr << "  Error opening file "<< fichentra2 <<"\n";
            return 1;}
        nlin2=0;
        while (getline(entrada, line)){
            std::istringstream iss(line);
            if (!(iss >> c >> d2o >> d2p)){
                std::cerr << " Format error in file "<< fichentra2<<std::endl;
                entrada.close();
                return 1;}
            crec[nlin2]=c;
            d2obs[nlin2]=d2o;
            d2prd[nlin2]*=pow(d2p,double(1.0/repeticiones));
            ++nlin2;
            if (nlin2>=maxlin){break;}
        }
        entrada.close();
    }
    for (i=0;i<nlin;++i){
        NeH[i]=repeticiones/NeH[i];
        NeA[i]=NeA[i]/repeticiones;}
    for (i=0;i<nlin2;++i){
        dif2+=pow((d2obs[i]-d2prd[i]),2);}
    dif2/=nlin2;
    //Ne:
    std::string fichsal=fichero+"_Ne_estimates";
    salida.open(fichsal, std::ios::out);
    salida << "Ne averages over "<<repeticiones<<" independent estimates."<<"\n";
//    salida << "Generation"<<"\t"<<"Arithmetic_mean"<<"\t"<<"Geometric_mean"<<"\t"<<"Harmonic_mean"<<"\n";
    salida << "Generation"<<"\t"<<"Geometric_mean"<<"\n";
    for (i=0;i<nlin;++i){
//        salida << generacion[i]<<"\t"<<NeA[i]<<"\t"<<NeG[i]<<"\t"<<NeH[i]<<"\n";
        salida << generacion[i]<<"\t"<<NeG[i]<<"\n";
    }
    salida.close();
    //d2:
    std::string fichsal2=fichero+"_d2_sample";
    salida.open(fichsal2, std::ios::out);
    salida << "Sample d2 values. Average (d2obs-d2prd)^2 = " <<dif2<<"\n";
    salida << "rec.rate_c"<<"\t"<<"observed_d2"<<"\t"<<"adjusted_d2"<<"\n";
    for (i=0;i<nlin2;++i){
        salida << crec[i]<<"\t"<<d2obs[i]<<"\t"<<d2prd[i]<<"\n";
    }
    salida.close();

    return 0;
}

