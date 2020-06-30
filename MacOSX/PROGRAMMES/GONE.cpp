//
//  main.cpp
//  GONE (Genetic Optimization for Ne Estimation)
//
//  Created by Enrique Santiago on 12/11/19.
//  Copyright © 2019 Enrique Santiago. All rights reserved.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <ctime>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <random>
#include <limits>

#define nlinmax 500
#define ngenmax 2000

void CalculaSC();

double SC;
double nbin[nlinmax];
double cval[nlinmax];
double d2cobs[nlinmax];
double xnbin[nlinmax];
double xcval[nlinmax];
double xd2cobs[nlinmax];
int indx[nlinmax];
double cval12[nlinmax];
double cvalrep[nlinmax];
double cval1masc2[nlinmax];
double cval2[nlinmax];
double d2cprd[nlinmax];
double Neblock[nlinmax];
int segdesdehasta[nlinmax];
double Ne[ngenmax];
double sumNe[ngenmax];
double d2c;
double Necons, Nemed;
int nlin, gmax,gmax2,gmax3, nsegmentos;
int muestrasalida=10;
int hp=2;
double fval,fvalsample, correccion=1;
double samplex,sampley,sampleZ2,sampleZ3,sampleZ4;
int conta,conta2;
double sample_size=0, sample_size_h=0;
bool flagporbloques=false, flagsolape=true;
double maxdouble=std::numeric_limits< double >::max();
unsigned long long int semilla=0;
int resolucion, topeposigen,ndes,gen, maxgen,des, ind1,ind2, ind3, posigen, posiblock, ancho1, ancho2;
int ndesini, indexmaxSC,indexminSC, maxsegmentos, posi1,posi2,posi3, contagen, topeposiinicio, topeposiinicio1,topeposiinicio4;
int tercio1, tercio2, tercio12, nhijos,topeposiinicio2,topeposiinicio3;
double efectomut,efectomutlateral, frecmut, frecrec,frecnorec, frecinversion, frecmutlateral, efecto, maxSC, minSC;
double efectomutsuave, SCmed, SCmed2, nsegmed,nsegmed2, SCbest,nsegbest, SCmedanterior,increfval,frecmutfval;
double tope,topesalto,invtopesalto,topesalto2,invtopesalto2;
bool flag;
struct serie { //++
    double efval; //++
    double SCval; //++
    int nseg; //++
    int segbl[nlinmax]; //++
    double Nebl[nlinmax]; //++
}; //++
serie bichoP[1000];//++ Parents
serie bichoH[1000];//++ Offspring
serie bicho;//++

std::random_device seed;  //used to obtain a seed for the random number engine
std::mt19937_64 genera(seed()); // mersenne_twister_engine 64bit (very good but very big)
std::uniform_real_distribution<> uniforme01(0.0, 1.0); //uses the result of the engine to generate uniform dist

int main(int argc, char* argv[]){
    double a,aa,b,bb,d,maxc, Ne0, increm, producto=1,mincval;
    double clow=0, chigh=0.5, minn, sumnbins=0,dpob;
    int i,ii,jj, j,imaxc;
    int mini, sizebins=1,nbins=50;
    bool  resizebins=false, flagdebug=false, flagrep=false,flagne=false;
    bool lcflag=false,hcflag=false, flagwrong=false, flagnbins=false;
    struct tm *loctime;
    long tinicio, tfinal;
    clock_t tini, tiniabs, tpas;
    
    if ((argc<2)){
        std::cerr << " Usage: "<< argv[0]<< " [-lc lowest_c_value] [-hc highest_c_value] [-ne number_of_estimates] [-sd seed] [-sr] [-vb] input_file_name [output_file_name]"<<std::endl;
        return 1;}
    
    for (i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            std::cout <<std::endl;
            std::cout << "  GONE v1.0. Genetic Optimization for Ne Estimation."<<std::endl<<std::endl;
            std::cout << "  Usage: "<< argv[0]<< " [-lc lowest_c_value] [-hc highest_c_value] [-ne number_of_estimates] [-sd seed] [-sr] [-vb] input_file_name [output_file_name]"<<std::endl<<std::endl;
            std::cout << "  The program searches for the best series of historical Ne values for a"<<std::endl;
            std::cout << "    series of observed LD values between pairs of loci. The first line in"<<std::endl;
            std::cout << "    the input file is a code (2, 1 or 0) which refers to the genetic data"<<std::endl;
            std::cout << "    type (diploid unphased, phased or pseudohaploid). The second line in"<<std::endl;
            std::cout << "    the input file contains the sample size. The third line is the devia-"<<std::endl;
            std::cout << "    tion from Hardy-Weimberg equilibrium, that is the correlation f of"<<std::endl;
            std::cout << "    alleles within locus. Each of the following lines corresponds to a"<<std::endl;
            std::cout << "    particular bin of recombination (c) and contains the number of SNP"<<std::endl;
            std::cout << "    pairs in the bin, the harmonic mean of c values and the arithemetic "<<std::endl;
            std::cout << "    mean ofthe observed d2 values separated by tabs or spaces. Additional"<<std::endl;
            std::cout << "    columns are ignored. Comment lines are preceded by a # symbol."<<std::endl;
            std::cout << "  Switches:"<<std::endl;
            std::cout << "   -lc  lowest recombination frequency to be considered."<<std::endl;
            std::cout << "        If not given, lc is set to 0.001."<<std::endl;
            std::cout << "   -hc  highest recombination frequency to be considered."<<std::endl;
            std::cout << "        If not given, hc is set to 0.5."<<std::endl;
            std::cout << "   -ne  number of the best Ne estimates to be averaged (default 10)"<<std::endl;
            std::cout << "   -sd  integer to seed random number generator. Taken "<<std::endl;
            std::cout << "        from the system, if not given."<<std::endl;
            std::cout << "   -sr  samplig with replacement."<<std::endl;
            std::cout << "   -vb  verbose mode with several intermediate results."<<std::endl;
            std::cout << "  Four files are generated with the following extensions:"<<std::endl;
            std::cout << "   _GONE_Ne Estimates of Ne (col 2) backward in time."<<std::endl;
            std::cout << "   _GONE_d2 Observed d2 (col 2), expected d2 in the sample(col 3)"<<std::endl;
            std::cout << "    and c values (col 1)."<<std::endl;
            std::cout << "   _GONE_input The input file rearranged."<<std::endl;
            std::cout << "   _GONE_log Log file."<<std::endl;
            return 0;}
        if (arg == "-vb"|| arg== "--verbose") {
            flagdebug=true;}
        if (arg == "-bs"|| arg== "--bin-size") {
            resizebins=true;
            sizebins=std::atoi(argv[i+1]);
            if (sizebins<1000){
                std::cerr <<" Invalid number of pairs of SNPs. It must be larger than 1000."<< std::endl;
                return 1;}}
        if (arg == "-ne"|| arg== "--number-of-estimates") {
            flagne=true;
            muestrasalida=std::atoi(argv[i+1]);
            if (muestrasalida<1 || muestrasalida>20){
                std::cerr <<" Invalid number of estimates to be averaged. It must be from 1 to 20."<< std::endl;
                return 1;}}
        if (arg == "-bn"|| arg== "--bin-number") {
            flagnbins=true;
            nbins=std::atoi(argv[i+1]);
            if (nbins<15){
                std::cerr <<" Invalid number of bins. It must be larger than 14."<< std::endl;
                return 1;}}
        if (arg == "-sr") {
            flagrep=true;}
        if (arg == "-lc"|| arg== "--lowest-c") {
            lcflag=true;
            clow=std::atof(argv[i+1]);
            if (clow<0 || clow>0.5){
                std::cerr <<" Invalid recombination frequency. It must be between 0 and 0.5."<< std::endl;
                return 1;}}
        if (arg == "-hc"|| arg== "--highest-c") {
            hcflag=true;
            chigh=std::atof(argv[i+1]);
            if (chigh<0 || chigh>0.5){
                std::cerr <<" Invalid recombination frequency. It must be a number between 0 and 0.5."<< std::endl;
                return 1;}}
        if (arg == "-sd"|| arg== "--seed") {
            semilla=std::atoi(argv[i+1]);}
    }
    if (chigh<=clow){
        std::cerr <<" Invalid range of recombination frequencies."<< std::endl;
        return 1;}
    
    // ***************************************************
    i=1;
    if (resizebins){i=i+2;}
    if (flagne){i=i+2;}
    if (flagnbins){i=i+2;}
    if (flagrep){i=i+1;}
    if (hcflag){i=i+2;}
    if (lcflag){i=i+2;}
    if (semilla!=0){i=i+2;
        genera.seed(semilla);}
    if (flagdebug){i=i+1;}
    std::string fichero=argv[i];
    std::string fichsal=argv[i];
    if (argc==(i+2)){
        fichsal=argv[i+1];}
    
    //            std::string fichero="ISLAS_KANK2";
    //            std::string fichsal="ISLAS_KANK2";
    
    std::string fichero_sal_NeH=fichsal+"_GONE_Nebest";
    std::string fichero_sal_d2=fichsal+"_GONE_d2";
    std::string fichero_sal_log=fichsal+"_GONE_log";
    std::string fichero_sal_rearrange=fichsal+"_GONE_input";
    std::string fichero_sal_evol=fichsal+"_GONE_evol";   // --debug
    std::string fichero_caca=fichsal+"_GONE_caca";   // --debug
    std::ofstream salida;
    
    if (!lcflag){clow=0.001;} //
    tinicio=time(NULL);
    loctime=localtime(&tinicio);
    std::ifstream entrada;
    entrada.open(fichero, std::ios::in);
    if (!(entrada.is_open())){
        std::cerr << " Error opening file "<< fichero <<"\n";
        return 1;}
    nlin=-3;
    std::string line;
    while (getline(entrada, line)){
        if (line[0]!='#' && line[0]!='*' && line[0]!='>' && line[0]!='/'){
            std::istringstream iss(line);
            if (nlin==-3){
                if (!(iss >> a)){
                    std::cerr << " Format error in file "<< fichero<<std::endl;
                    entrada.close();
                    return 1;
                    break; }
                if (a!=0 && a!=1 && a!=2){
                    std::cerr <<" Specify a valid code for type of genotyping data (0, 1 or 2)."<< std::endl;
                    return 1;}
                hp=(int)(a);
                nlin++;
            }
            else if (nlin==-2){
                if (!(iss >> a)){
                    std::cerr << " Format error in file "<< fichero<<std::endl;
                    entrada.close();
                    return 1;
                    break; }
                if (a<2){
                    std::cerr <<" Specify a sample size larger than 1."<< std::endl;
                    return 1;}
                sample_size=a;
                nlin++;
            }
            else if (nlin==-1){
                if (!(iss >> a)){
                    std::cerr << " Format error in file "<< fichero<<std::endl;
                    entrada.close();
                    return 1;
                    break; }
                if (a<-1 || a>1){
                    std::cerr <<" Specify a f value in the range -1 and +1."<< std::endl;
                    return 1;}
                fvalsample=a;
                nlin++;
            }
            else{
                if (!(iss >> a >> b >> d)) {
                    std::cerr << " Format error in file "<< fichero<<std::endl;
                    entrada.close();
                    return 1;
                    break; } // error
                if ((a<0) || (b<=0) || (b>0.5)){flagwrong=true;}
                if ((b>=clow) && (b<=chigh) && (a>0)){
                    xnbin[nlin]=a;
                    xcval[nlin]=b;
                    xd2cobs[nlin]=d;
                    indx[nlin]=nlin;
                    nlin++;}
            }
            if (nlin==nlinmax){
                break;}
        }
    }
    
    for (i=0;i<nlin;++i){
        d2cobs[i]=0;
        cval[i]=0;
        nbin[i]=0;
    }
    
    sample_size_h=sample_size*2;
    samplex=(double)((sample_size_h-2.0)*(sample_size_h-2.0)*(sample_size_h-2.0)+8.0/5.0*(sample_size_h-2.0)*(sample_size_h-2.0)+4*(sample_size_h-2.0))/((sample_size_h-1.0)*(sample_size_h-1.0)*(sample_size_h-1.0)+(sample_size_h-1.0)*(sample_size_h-1.0));
    sampley=(double)(2.0*sample_size_h-4.0)/((sample_size_h-1.0)*(sample_size_h-1.0));
    sampleZ3=(double)(1.0)/(sample_size-1.0);
    sampleZ4=(double)(1.0)/(sample_size_h-1.0);
    sampleZ2=(double)1.0-0.8/(sample_size-1.0);
    if (flagwrong){
        std::cerr <<" Wrong data type. Probably some input values are out of range."<< std::endl;
        return 1;}
    // Sort the input by c values and group them if xnbin<sizebins
    for(i=0;i<(nlin-1);++i){
        imaxc=i;
        maxc=xcval[indx[imaxc]];
        for (j=i+1;j<nlin;++j){
            if (maxc<xcval[indx[j]]){
                imaxc=j;
                maxc=xcval[indx[imaxc]];
            }}
        jj=indx[i];
        indx[i]=indx[imaxc];
        indx[imaxc]=jj;
    }
    
    //MUESTREA LAS PRIMERAS LINEAS PARA HACER UN TANTEO DE Ne
    Nemed=0;
    conta=0;
    for (i=5;i<25;++i){
        if (i>=nlin){break;}
        producto=(1.0-xcval[indx[i]])*(1-xcval[indx[i]]);
        Ne0=1/((xd2cobs[indx[i]]-1/sample_size)/producto*(1-producto))-2.2*producto/(1-producto);
        if (Ne0>10){
            conta+=1;
            Nemed+=(Ne0);
        }
    }
    if (conta>0){Nemed=(Nemed/conta);}
    else{Nemed=1000;}
    if (Nemed<10){Nemed=10;}
    //MIRA LA GENERACION MAXIMA QUE TODAVIA CONTIENE UN 5% DE INFORMACION
    gmax=int(log10(0.05)/log10(1.0-1.0/Nemed)+.5);
    mincval=clow;
    
    // VA REDUCIENDO HASTA QUE NO HAYA LINEAS CON MENOS DE UN No DE PAREJAS DE SNPs
    if (!resizebins){sizebins=10000;}
    while (nlin>18){
        // Busca el mínimo
        mini=0;
        minn=xnbin[indx[mini]];
        for (i=1;i<nlin;++i){
            if (minn>xnbin[indx[i]]){mini=i;minn=xnbin[indx[i]];}}
        if (minn>sizebins){break;}
        // Mira cual de los adyacentes es el menor
        j=mini-1;
        jj=mini+1;
        if (j<0){j=jj;}
        else{
            if (jj>=nlin){j=j;}
            else{
                if (xnbin[indx[j]]>xnbin[indx[jj]]){j=jj;}}}
        if (j<mini){mini=j;j+=1;}
        //los suma
        increm=(xnbin[indx[mini]]/xcval[indx[mini]]+xnbin[indx[j]]/xcval[indx[j]])/(xnbin[indx[mini]]+xnbin[indx[j]]);
        xcval[indx[mini]]=1.0/increm;
        increm=(xnbin[indx[mini]]*xd2cobs[indx[mini]]+xnbin[indx[j]]*xd2cobs[indx[j]])/(xnbin[indx[mini]]+xnbin[indx[j]]);
        xd2cobs[indx[mini]]=increm;
        xnbin[indx[mini]]=xnbin[indx[mini]]+xnbin[indx[j]];
        //comprime la lista
        --nlin;
        for (i=j;i<nlin;++i){
            xcval[indx[i]]=xcval[indx[i+1]];
            xd2cobs[indx[i]]=xd2cobs[indx[i+1]];
            xnbin[indx[i]]=xnbin[indx[i+1]];}
    }
    if (!flagnbins){
        // CALCULA EL No DE LINEAS FINALES A PARTIR DEL NUMERO INICIAL
        producto=log10(xnbin[indx[nlin-1]]);
        producto=producto-3.0;
        if (producto<0){producto=0;}
        producto=pow(2.0,producto)*10;
        nbins=int(producto+8);
        if((producto>nlin)){nbins=nlin;}
        if((producto<30)){nbins=30;}
        if (nbins>60){nbins=60;}}  //LIMITE EN EL NUMERO DE LINEAS
    
    //VA REDUCIENDO HASTA DEJAR EL No DE LINEAS
    while(nlin>nbins){
        // Busca el mínimo
        mini=0;
        minn=xnbin[indx[mini]];
        for (i=1;i<nlin;++i){
            if (minn>xnbin[indx[i]]){mini=i;minn=xnbin[indx[i]];}}
        // Mira cual de los adyacentes es el menor
        j=mini-1;
        jj=mini+1;
        if (j<0){j=jj;}
        else{
            if (jj>=nlin){j=j;}
            else{
                if (xnbin[indx[j]]>xnbin[indx[jj]]){j=jj;}}}
        if (j<mini){mini=j;j+=1;}
        //los suma
        increm=(xnbin[indx[mini]]/xcval[indx[mini]]+xnbin[indx[j]]/xcval[indx[j]])/(xnbin[indx[mini]]+xnbin[indx[j]]);
        xcval[indx[mini]]=1.0/increm;
        increm=(xnbin[indx[mini]]*xd2cobs[indx[mini]]+xnbin[indx[j]]*xd2cobs[indx[j]])/(xnbin[indx[mini]]+xnbin[indx[j]]);
        xd2cobs[indx[mini]]=increm;
        xnbin[indx[mini]]=xnbin[indx[mini]]+xnbin[indx[j]];
        //comprime la lista
        --nlin;
        for (i=j;i<nlin;++i){
            xcval[indx[i]]=xcval[indx[i+1]];
            xd2cobs[indx[i]]=xd2cobs[indx[i+1]];
            xnbin[indx[i]]=xnbin[indx[i+1]];}
    }
    
    if (nlin<15){
        std::cerr <<" Too few bins."<< std::endl;
        return 1;}
    
    mincval=maxdouble;
    fval=(1.0+fvalsample*(sample_size_h-1.0))/(sample_size_h-1+fvalsample);
    if (hp==2){correccion=1.0/((1.0+fval)*(1.0+fval));}
    else if (hp==1){correccion=1;}
    else {correccion=1;}
    sumnbins=0;
    for (i=0;i<nlin;++i){
        nbin[i]=xnbin[indx[i]];
        sumnbins+=nbin[i];
        cval[i]=xcval[indx[i]];
        
        if (cval[i]<mincval){mincval=cval[i];}
        d2cobs[i]=xd2cobs[indx[i]];
        cval12[i]=(1-cval[i])*(1-cval[i]);
        if (flagrep){
            cvalrep[i]=1.0;}
        else{
            cvalrep[i]=cval12[i];}
        cval1masc2[i]=(1+cval[i]*cval[i]);
        cval2[i]=(cval[i]*cval[i]);
    }
    conta=0;
    Nemed=0;
    for (i=2;i<8;++i){
        if (d2cobs[i]>0){
            if (hp==0){dpob=4.0*(d2cobs[i]-sampleZ3)/(cvalrep[i]*sampleZ2);}
            else if (hp==1){dpob=(d2cobs[i]-sampleZ4)/(cvalrep[i]);}
            else {dpob=(d2cobs[i]/((1+fval)*(1+fval))-sampley)/(samplex*cvalrep[i]);}
            Ne0=2*(cval1masc2[i]/(dpob*2*(1-cval12[i]))-1.1*cval12[i]/(1-cval12[i]));
            if (Ne0>10){
                Nemed=Nemed+Ne0;
                conta=conta+1;}}
    }
    if (conta>0){Nemed=(Nemed/conta);}
    else {Nemed=2000;}
    if (Nemed<10){Nemed=10;}
    entrada.close();
    tiniabs=clock();
    tpas=clock() - tiniabs;
    salida.open(fichero_sal_log, std::ios::out);
    salida << "(GONE v1.0)\n\n";
    salida << "Input File: "<< fichero <<"\n\n";
    salida << "Command:";
    for (i=0;i<argc;++i){salida << " " << argv[i];}
    salida << "\n\n";
    salida << "START: " << asctime(loctime)<< "\n";
    salida.close();
   
    if (flagdebug){
        salida.open(fichero_sal_evol, std::ios::out);
        salida << "Gener"<<"\t"<<"SCbest"<<"\t"<<"SCmed1"<<"\t"<<"nsegbest"<<"\t"<<"nsegmed"<<"\n";
        salida.close();}
    
    tinicio=time(NULL);
    tini=clock();
    std::cout << "  Start of processing: "<<fichsal<< std::endl;
    
    gmax=int((1.0/(mincval)));
    if (gmax>(ngenmax-10)){gmax=ngenmax-10;}
    gmax2=int(gmax*26/40);
    gmax3=gmax2+int((gmax-gmax2)/2);
    
    //                                    PARAMETROS DE AJUSTE
    //**************************************************************************************
    //**************************************************************************************
    //**************************************************************************************
    flagsolape=true;
    flagporbloques=true;
    resolucion=3;
    topeposigen=gmax2-2*resolucion;
    topeposiinicio1=60;
    topeposiinicio2=120;
    topeposiinicio3=240;
    topeposiinicio4=480;
    if (topeposiinicio1>topeposigen){topeposiinicio1=topeposigen;}
    if (topeposiinicio2>topeposigen){topeposiinicio2=topeposigen;}
    if (topeposiinicio3>topeposigen){topeposiinicio3=topeposigen;}
    if (topeposiinicio4>topeposigen){topeposiinicio4=topeposigen;}
    if (topeposiinicio2>topeposigen){topeposiinicio2=topeposigen;}
    topeposiinicio=topeposigen;
    topesalto=9.0;
    invtopesalto=1.0/topesalto;
    topesalto2=50;
    invtopesalto2=1.0/topesalto2;
    maxgen=750;  //
    ndesini=3000;
    ndes=1000;
    increfval=0.001;
    frecmutfval=0.0;
    efectomutsuave=0.02;
    frecinversion=0.3; //INVERSION DE LOS DOS ULTIMOS
    tercio1=10; // aqui van los mejores reproductores que no se tocan
    tercio2=90; // son los siguientes y se sustituyen por hijos si estos son mejores
    tercio12=tercio1+tercio2;
    nhijos=tercio12;
    
    
    // PRECARGA ALEATORIA CON SELECCION
    bicho.efval=fval; //++
    bicho.nseg=4; //++
    bicho.segbl[0]=0; //++
    bicho.segbl[3]=gmax3; //++
    bicho.segbl[4]=gmax; //++
    for (ii=0;ii<ndesini;++ii){
        for (;;){
            aa=uniforme01(genera);
            if (aa<0.4){
                posi1=int(uniforme01(genera)*(topeposiinicio2-resolucion))+resolucion;
                posi2=int(uniforme01(genera)*(topeposiinicio2-resolucion))+resolucion;}
            else if (aa<0.8){
                posi1=int(uniforme01(genera)*(topeposiinicio3-resolucion))+resolucion;
                posi2=int(uniforme01(genera)*(topeposiinicio3-resolucion))+resolucion;}
            else{
                posi1=int(uniforme01(genera)*(topeposiinicio4-resolucion))+resolucion;
                posi2=int(uniforme01(genera)*(topeposiinicio4-resolucion))+resolucion;}
            
            if (posi1>posi2){posi3=posi1;posi1=posi2;posi2=posi3;}
            if ((posi2-posi1)>5){break;}}
        bicho.segbl[1]=posi1; //++
        bicho.segbl[2]=posi2; //++
        bicho.Nebl[0]=Nemed;
        tope=uniforme01(genera)*2*topesalto;
//        tope=topesalto;  // ALTERNATIVA A LA LINEA ANTERIOR
        for (jj=1;jj<3;++jj){
            aa=1+uniforme01(genera)*(tope); //ERA2
            if (uniforme01(genera)<0.66){aa=1/aa;} //ERA0.5
            aa=Nemed*aa;
            if (aa>10000000){aa=10000000;}
            if (aa<5){aa=5;}
            bicho.Nebl[jj]=aa; //++
        }
        aa=1.0+uniforme01(genera)*tope/2; //era 0.5
        aa=Nemed*aa;
        if (aa>10000000){aa=10000000;}
        if (aa<5){aa=5;}
        bicho.Nebl[3]=aa; //++
        CalculaSC();
        bicho.SCval=SC; //++
        // Selecciona los padres
        if (ii<tercio12){
            bichoP[ii]=bicho;}
        else{
            maxSC=0; //Busca el peor SC entre los futuros padres
            indexmaxSC=0;
            for (jj=0;jj<tercio12;++jj){
                if (bichoP[jj].SCval>maxSC){indexmaxSC=jj;maxSC=bichoP[jj].SCval;}
            }
            if (bicho.SCval<bichoP[indexmaxSC].SCval){
                bichoP[indexmaxSC]=bicho;}
        }
    }
    // Ordena los padres resultantes
    for (ii=0;ii<tercio12-1;++ii){
        indexminSC=ii;
        for (jj=ii+1;jj<tercio12;++jj){
            if (bichoP[jj].SCval<bichoP[indexminSC].SCval){indexminSC=jj;}}
        if (indexminSC!=ii){//swap
            bicho=bichoP[ii];
            bichoP[ii]=bichoP[indexminSC];
            bichoP[indexminSC]=bicho;
        }
    }
    //              COMIENZAN LOS CICLOS
    contagen=0;
    SCmedanterior=SC;
    if(flagdebug){
        salida.open(fichero_caca, std::ios::app);
        for (j=0;j<bichoP[0].nseg;++j){
            salida << j<<"\t"<<bichoP[0].segbl[j] <<"\t" << bichoP[0].Nebl[j]/2.0 <<"\n";}
        salida<<"\n";
        salida.close();}
    for (gen=0;gen<maxgen;++gen){
        
        if(flagdebug){
            if (gen==int(double(gen/100.0))*100){
                salida.open(fichero_caca, std::ios::app);
                for (j=0;j<bichoP[0].nseg;++j){
                    salida <<"-> "<< gen<<"\t"<< j<<"\t"<<bichoP[0].segbl[j] <<"\t" << bichoP[0].Nebl[j]/2.0 <<"\n";}
                salida<<"\n";
                salida.close();}}
        
        switch (gen){
            case 0:{frecmut=0.3;efectomut=0.3;efectomutlateral=5;frecrec=0.6;frecnorec=0.3;
                frecmutlateral=0.2;maxsegmentos=4;flagsolape=true;break;}
            case 300:{frecmut=0.2;efectomut=0.2;efectomutlateral=2;frecrec=0.2;frecnorec=0.5;
                frecmutlateral=0.2;maxsegmentos+=2;flagsolape=true;break;}
            case 600:{frecmut=0.5;efectomut=0.05;efectomutlateral=1;frecrec=0.2;frecnorec=0.5;
                frecmutlateral=0.2;maxsegmentos+=10;flagsolape=true;break;}
            case 700:{frecmut=0.5;efectomut=0.2;efectomutlateral=1;frecrec=1;frecnorec=0;
                frecmutlateral=0.2;maxsegmentos+=10;flagsolape=false;break;}
            case 710:{frecmut=0.5;efectomut=0.2;efectomutlateral=1;frecrec=1;frecnorec=0;
                frecmutlateral=0.2;maxsegmentos+=10;flagsolape=false;break;}
            case 720:{frecmut=0.5;efectomut=0.1;efectomutlateral=1;frecrec=0.95;frecnorec=0;
                frecmutlateral=0.2;maxsegmentos+=10;flagsolape=false;break;}
            case 730:{frecmut=0.5;efectomut=0.04;efectomutlateral=1;frecrec=0.95;frecnorec=0;
                frecmutlateral=0.2;maxsegmentos+=20;flagsolape=false;break;}
        }

        for (des=0;des<ndes;++des){
            for (;;){
                for (;;){
                    flag=true;
                    ind1=int(uniforme01(genera)*tercio1);
                    ind2=int(uniforme01(genera)*tercio12);
                    if (uniforme01(genera)<0.5){ind3=ind1;ind1=ind2;ind2=ind3;}
                    // Recombinacion ENTRE y DENTRO CON IGUAL PROBABILIDAD: genera una serie a partir de dos o un padre
                    if (uniforme01(genera)<frecrec){
                        aa=uniforme01(genera);
                        if (aa<0.1) {topeposiinicio=topeposiinicio1;}
                        else if (aa<0.5) {topeposiinicio=topeposiinicio2;}
                        else if (aa<0.9) {topeposiinicio=topeposiinicio3;}
                        else {topeposiinicio=topeposiinicio4;}
                        if (flagporbloques){
                            do{
                                posiblock=int(uniforme01(genera)*(bichoP[ind1].nseg-1));
                                posigen= bichoP[ind1].segbl[posiblock]+int(uniforme01(genera)*(bichoP[ind1].segbl[posiblock+1]-bichoP[ind1].segbl[posiblock]));
                            }
                            while (posigen>topeposiinicio);
                        }
                        else{
                            posigen=resolucion+int(uniforme01(genera)*topeposiinicio);
                        }
                        if (uniforme01(genera)<.5){  //DENTRO: iguala individuos y centra posigen
                            ind2=ind1;
                        }
                        for (jj=0;jj<=bichoP[ind1].nseg;++jj){
                            if (abs(bichoP[ind1].segbl[jj]-posigen)<=resolucion){flag=false;break;}
                        }
                        if (flag){
                            for (jj=0;jj<=bichoP[ind2].nseg;++jj){
                                if (abs(bichoP[ind2].segbl[jj]-posigen)<=resolucion){flag=false;break;}
                            }
                        }
                        if (flag){
                            ii=0;
                            for (jj=0;jj<bichoP[ind1].nseg;++jj){
                                if (bichoP[ind1].segbl[jj]<posigen){
                                    bicho.segbl[ii]=bichoP[ind1].segbl[jj];
                                    bicho.Nebl[ii]=bichoP[ind1].Nebl[jj];
                                    ++ii;}
                                else{
                                    break;
                                }
                            }
                            for (jj=1;jj<bichoP[ind2].nseg+1;++jj){
                                if (bichoP[ind2].segbl[jj]>posigen){
                                    bicho.segbl[ii]=posigen;
                                    bicho.Nebl[ii]=bichoP[ind2].Nebl[jj-1];
                                    ++ii;
                                    break;
                                }
                            }
                            for (jj=1;jj<bichoP[ind2].nseg;++jj){
                                if (bichoP[ind2].segbl[jj]>posigen){
                                    bicho.segbl[ii]=bichoP[ind2].segbl[jj];
                                    bicho.Nebl[ii]=bichoP[ind2].Nebl[jj];
                                    ++ii;
                                }
                            }
                            bicho.nseg=ii;
                            bicho.segbl[ii]=bichoP[ind2].segbl[bichoP[ind2].nseg];
                        }
                    }
                    else{bicho=bichoP[ind1];}
                    if (flag){break;}
                }
                // Fusion de dos bloques con probabilidad frecnorec haciendo la media harmonica de los Ne
                if (bicho.nseg>3){
                    if ((uniforme01(genera)<frecnorec) || (bicho.nseg>maxsegmentos)){
                        posiblock=int(uniforme01(genera)*(bicho.nseg-2)); //FUSIONA AL AZAR
                        ancho1=bicho.segbl[posiblock+1]-bicho.segbl[posiblock];
                        ancho2=bicho.segbl[posiblock+2]-bicho.segbl[posiblock+1];  // Media harmonica
                        aa=(ancho1+ancho2)/(ancho1/bicho.Nebl[posiblock]+ancho2/bicho.Nebl[posiblock+1]);
                        flag=true;
                        if (posiblock>0){
                            bb=aa/bicho.Nebl[posiblock-1];
                            //                            if (bb<invtopesalto || bb>topesalto){flag=false;}}
                            if (bb<invtopesalto2){flag=false;}}
                        if (posiblock<bicho.nseg-2){
                            bb=aa/bicho.Nebl[posiblock+2];
                            //                            if (bb<invtopesalto || bb>topesalto){flag=false;}}
                            if (bb>topesalto2){flag=false;}}
                        if (flag){
                            bicho.Nebl[posiblock]=aa;
                            for (jj=posiblock+1;jj<=bicho.nseg;++jj){
                                bicho.Nebl[jj]=bicho.Nebl[jj+1];
                                bicho.segbl[jj]=bicho.segbl[jj+1];
                            }
                            --bicho.nseg;
                        }
                    }
                }
                // Control del tope de número de segmentos
                flag=true;
                if ((bicho.nseg>maxsegmentos) || (bicho.nseg<3)){flag=false;}
                if (flag){break;}
            }
            // Mutación: Cambia el valor de ne de segmentos al azar
            for (posiblock=0;posiblock<bicho.nseg;++posiblock){
                if (uniforme01(genera)<frecmut){
                    efecto= bicho.Nebl[posiblock]*uniforme01(genera)*efectomut;}
                else{
                    efecto= bicho.Nebl[posiblock]*uniforme01(genera)*efectomutsuave;}
                if (uniforme01(genera)<0.5){efecto=-efecto;}
                aa=bicho.Nebl[posiblock]+efecto;
                flag=true;
                if (posiblock>0){
                    bb=aa/bicho.Nebl[posiblock-1];
                    if (bb<invtopesalto2){flag=false;}}
                if (posiblock<bicho.nseg-1){
                    bb=aa/bicho.Nebl[posiblock+1];
                    if (bb>topesalto2){flag=false;}}
                if (flag){
                    if (aa>10000000){aa=10000000;}
                    if (aa<5){aa=5;}
                    bicho.Nebl[posiblock]=aa;}
            }
            // Iguala dos dos ultimos:
            aa=bicho.Nebl[bicho.nseg-1];
            bb=bicho.Nebl[bicho.nseg-2];
            if (aa>(bb*1.2)){bicho.Nebl[bicho.nseg-1]=bb*1.2;}
            if (aa<(bb/1.2)){bicho.Nebl[bicho.nseg-1]=bb/1.2;}
            
            for (posiblock=0;posiblock<bicho.nseg;++posiblock){
                if (bicho.Nebl[posiblock]<5){bicho.Nebl[posiblock]=5;}}
            
            // Mutación lateral: Cambia el límite de segmentos al azar
            for (posiblock=1;posiblock<bicho.nseg-1;++posiblock){
                if (uniforme01(genera)<frecmutlateral){
                    efecto= int(uniforme01(genera)*efectomutlateral)+1;}
                else {
                    efecto=1;
                    if (uniforme01(genera)<0.5){efecto=0;}}
                if (uniforme01(genera)<0.5){ //ERA0.6
                    if ((bicho.segbl[posiblock]-bicho.segbl[posiblock-1])>resolucion+efecto){
                        bicho.segbl[posiblock]-=efecto;}}
                else{
                    if ((bicho.segbl[posiblock+1]-bicho.segbl[posiblock])>resolucion+efecto){
                        bicho.segbl[posiblock]+=efecto;}}
            }
            
            bicho.efval=bichoP[ind1].efval;
            CalculaSC();
            bicho.SCval=SC; //++
            if (des<nhijos){
                bichoH[des]=bicho;}
            else{
                maxSC=-1;
                indexmaxSC=0;
                for (i=0;i<nhijos;++i){
                    if (bichoH[i].SCval>maxSC){
                        maxSC=bichoH[i].SCval;
                        indexmaxSC=i;}}
                if (SC<maxSC){
                    bichoH[indexmaxSC]=bicho;}
            }
        }
        //El primer tercio de los bichoP no se toca.
        
        if (flagsolape){
            //Mete los mejores bichoH en el segundo tercio de los bichoP si es que son mejores
            for (ii=tercio1;ii<tercio12;++ii){
                flag=true;
                minSC=maxdouble; //Busca el mínimo SC entre los hijos
                indexminSC=0;
                for (i=0;i<nhijos;++i){
                    if (bichoH[i].SCval<minSC){
                        minSC=bichoH[i].SCval;
                        indexminSC=i;}}
                for (i=tercio1;i<tercio12;++i){ // mira si ese mínimo es mejor que algún bichoP del segundo tercio
                    if (minSC<bichoP[i].SCval){
                        bichoP[i]=bichoH[indexminSC];
                        bichoH[indexminSC].SCval=maxdouble;
                        flag=false;
                        break;
                    }
                }
                if (flag){break;}}}
        else{
            //COPIA LOS HIJOS EN LOS PADRES:
            for (ii=0;ii<tercio12;++ii){
                bichoP[ii]=bichoH[ii];}}
        
        // Inversión: Intercambia dos ultimosNe adyacentes con probabilidad frecinversion
        if (uniforme01(genera)<frecinversion){
            for (ii=0;ii<tercio12;++ii){
                bicho=bichoP[ii];
                if (bicho.nseg>=4){
                    posiblock=bicho.nseg-2;
                    aa=bicho.Nebl[posiblock+1];
                    flag=true;
                    bb=aa/bicho.Nebl[posiblock-1];
                    if (bb>invtopesalto && bb<topesalto){
                        bicho.Nebl[posiblock+1]=bicho.Nebl[posiblock];
                        bicho.Nebl[posiblock]=aa;
                        CalculaSC();
                        bicho.SCval=SC; //++
                        bichoP[ii]=bicho;
                    }
                }
            }
        }
        
        // Ordena los padres
        for (ii=0;ii<tercio12-1;++ii){
            indexminSC=ii;
            for (jj=ii+1;jj<tercio12;++jj){
                if (bichoP[jj].SCval<bichoP[indexminSC].SCval){indexminSC=jj;}}
            if (indexminSC!=ii){//swap
                bicho=bichoP[ii];
                bichoP[ii]=bichoP[indexminSC];
                bichoP[indexminSC]=bicho;
            }
        }
        
        SCmed=nsegmed=0;
        for (ii=0;ii<tercio12;++ii){
            SCmed+=bichoP[ii].SCval;
            nsegmed+=bichoP[ii].nseg;
        }
        SCmed/=tercio12;
        nsegmed/=tercio12;
        SCbest=bichoP[0].SCval;
        nsegbest=bichoP[0].nseg;
        if(flagdebug){
            salida.open(fichero_sal_evol, std::ios::app);
            salida << gen <<"\t" << SCbest<<"\t"<<SCmed<<"\t"<<nsegbest<<"\t"<<nsegmed <<"\n";
            salida.close();}
    }
    if(flagdebug){
        salida.open(fichero_caca, std::ios::app);
        for (j=0;j<bichoP[0].nseg;++j){
            salida << j<<"\t"<<bichoP[0].segbl[j] <<"\t" << bichoP[0].Nebl[j]/2.0 <<"\n";
        }
        salida.close();}
    
    bicho=bichoP[0]; // El mejor
    CalculaSC();//se cargan los Neblock, segdesdehasta,SC,...
    
    for (i=0;i<=bicho.nseg;++i){Neblock[i]=bicho.Nebl[i];segdesdehasta[i]=bicho.segbl[i];}
    
    conta=0;
    for (i=0;i<bicho.nseg;++i){
        for (j=segdesdehasta[i];j<segdesdehasta[i+1];++j){
            Ne[conta]=Neblock[i];
            conta=conta+1;}}
    conta2=conta;
    // MEDIA GEOMETRICA DE LOS muestrasalida PRIMEROS
    for (i=0;i<conta2;++i){
        sumNe[i]=1;}  //
    for (ii=0;ii<muestrasalida;++ii){
        conta=0;
        for (i=0;i<bichoP[ii].nseg;++i){
            for (j=bichoP[ii].segbl[i];j<bichoP[ii].segbl[i+1];++j){
                sumNe[conta]*=pow(bichoP[ii].Nebl[i],1.0/muestrasalida);
                conta=conta+1;
                if (conta>gmax){break;}
            }
            if (conta>gmax){break;}
        }
    }
    
    for (i=0;i<conta2;++i){
        Ne[i]=Ne[i]/2;sumNe[i]/=2;}  // Now, individuals are diploids (only for representation)
    
    tfinal=time(NULL);
    tpas=(clock() - tiniabs);
    salida.open(fichero_sal_NeH, std::ios::out);
    for (j=0;j<gmax3+1;++j){
        salida << j+1 << "\t" << sumNe[j]<< "\n";}
    salida.close();
    
    salida.open(fichero_sal_d2, std::ios::out);
    for (i=0;i<nlin;++i){
        salida << cval[i] << "\t" << d2cobs[i] << "\t" << d2cprd[i] << "\n";}
    salida.close();
    
    salida.open(fichero_sal_rearrange, std::ios::out);
    salida << hp<< "\n";
    salida << sample_size<< "\n";
    salida << fvalsample<< "\n";
    for (i=0;i<nlin;++i){
        salida  << nbin[i]<<"\t"<<cval[i]<<"\t"<< d2cobs[i]<< "\n";}
    salida.close();
    
    tfinal=time(NULL);
    tpas=(clock() - tiniabs);
    salida.open(fichero_sal_log, std::ios::app);
    salida << "END  : " << asctime(localtime(&tfinal)) <<"\n";
    salida << "Total run time in seconds: " << (float(tpas)/CLOCKS_PER_SEC) <<"\n\n";
    salida << "  Range of c values included in this analysis: " << mincval<<" to "<<chigh << "\n\n";
    salida << "  Sample size (diploids): " << sample_size <<"\n\n";
    salida << "  Average residual square (d2obs-d2prd)^2 per c bin: " << SC/double(nlin) << "\n\n";
    salida << "  Observed d2 values were transformed by the predicted factor 1/[1+f]^2= "<<correccion<<"\n\n";
    salida << "  Bins have been rearranged in groups and written to file :"<<fichero_sal_rearrange << "\n\n";
    salida << "  Estimates of Ne values for the best solution were written to : "<< fichero_sal_NeH << "\n\n";
    salida << "  Observed and predicted d2 values estimated over replications for each c bin were written to :"<<fichero_sal_d2 << "\n\n";
    salida.close();
    
    std::cout << "  End of processing: "<< fichsal <<"\n";
    
    return 0;}


// Core of computation:
// This subroutine calculates the values of SD2_c, SW_c and the predicted d2_c
// for a particular set of Ne_t values stored in Ne[].
void CalculaSC(){
    int i,ii, exponente,hastasegmento;
    double Ne12[nlinmax];
    double Ne122[nlinmax];
    double Ne121[nlinmax];
    double Ne102[nlinmax];
    double Nec122[nlinmax];
    double Nec121[nlinmax];
    double Nec102[nlinmax];
    double p1b,p1a, s1, r1a, SD2, SW, cval12aa, cval12ancho, Ne12ancho, Ne122ancho, Ne102ancho, A, B;
    double cval12base, Ne12base, Ne122base, Ne102base;
    double Aplus, Bplus;
    nsegmentos=bicho.nseg;
    for (i=0;i<=nsegmentos;++i){Neblock[i]=bicho.Nebl[i];segdesdehasta[i]=bicho.segbl[i];}
    Necons=Neblock[nsegmentos-1];
    SC=0;
    hastasegmento=nsegmentos;
    if (hastasegmento>1){--hastasegmento;}
    for (i=0;i<nsegmentos;++i){Ne121[i]=1-2.2/Neblock[i];Ne12[i]=1-2.0/Neblock[i];Ne122[i]=1-2.2/Neblock[i];Ne102[i]=1-0.2/Neblock[i];}
    for (ii=0;ii<nlin;++ii){
        for (i=0;i<nsegmentos;++i){Nec121[i]=cval12[ii]*Ne121[i];Nec102[i]=cval12[ii]*Ne102[i];Nec122[i]=cval12[ii]*Ne122[i];}
        p1b=p1a=1;
        s1=0;
        r1a=1;
        SD2=0;
        SW=0;
        cval12aa=1;
        for (i=0;i<hastasegmento;++i){
            exponente=segdesdehasta[i+1]-segdesdehasta[i];
            Ne12base=Ne12[i];
            Ne122base=Ne122[i];
            Ne102base=Ne102[i];
            cval12base=cval12[ii];
            Ne12ancho=Ne122ancho=Ne102ancho=cval12ancho=1;
            while (exponente>0){
                if (exponente & 1){
                    Ne12ancho*=Ne12base;
                    Ne122ancho*=Ne122base;
                    Ne102ancho*=Ne102base;
                    cval12ancho*=cval12base;}
                Ne12base*=Ne12base;
                Ne122base*=Ne122base;
                Ne102base*=Ne102base;
                cval12base*=cval12base;
                exponente>>=1;
            }
            A=(1-Ne122ancho*cval12ancho)/(1-(Nec122[i]));
            B=(1-Ne12ancho)/(1-Ne12[i]);
            SD2+=s1*p1a*B+p1b/Neblock[i]*cval12aa*(Ne12[i]*B-Nec122[i]*A)/(Ne12[i]-Nec122[i]);
            SW+=p1a*(1-Ne12ancho)/(2/Neblock[i]);
            s1+=r1a*cval12aa/Neblock[i]*(1-Ne102ancho*cval12ancho)/(1-Nec102[i]);
            r1a*=Ne102ancho;
            p1a*=Ne12ancho;
            p1b*=Ne122ancho;
            cval12aa*=cval12ancho;
        }
        Aplus=Nec122[nsegmentos-1]/(1-Nec122[nsegmentos-1]);
        Bplus=Ne12[nsegmentos-1]*Necons/2;
        SD2+=s1*p1a*Necons/2+p1b/Necons*cval12aa*(Bplus-Aplus)/(Ne12[nsegmentos-1]-Nec122[nsegmentos-1]);
        SW+=p1a*Necons/2;
        d2c=SD2/SW*(cval1masc2[ii]);  // Correccion acumulacion diploides *(1+c^2)
        
        if (hp==2){ // Correccion diploides unphased
            d2cprd[ii]=(d2c*cvalrep[ii]*samplex+cval2[ii]/Necons*samplex+sampley)/(correccion);}
        else if (hp==1){ // Correccion diploides phased
            d2cprd[ii]=(d2c*cvalrep[ii]+sampleZ4);}
        else if (hp==0){ // Correccion pseudohaploides
            d2cprd[ii]=(d2c/4*cvalrep[ii]*sampleZ2+sampleZ3);}
        
    }
    for (ii=0;ii<nlin;++ii){
        A=d2cobs[ii]-d2cprd[ii];
        SC+=(A*=A);
    }
}









