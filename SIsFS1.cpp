//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
////      _____    _    _____     _____   _____       ////
////      \  __\  |_|   \  __\   |  __/   \  __\      ////
////       \ \    | |    \ \     | |__     \ \        ////
////      __\ \   | |   __\ \    | __/    __\ \       ////
////      \____\  |_|   \____\   |_|      \____\      ////
////                                                  ////
//////////////////////////////////////////////////////////
////////////////////////by R3ZZ///////////////////////////
//////////////////////////////////////////////////////////

#include <cstring>
#include "stdafx.h"   
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include "SFS.h"
#include <cmath>
#include <iomanip>


#define pi 3.141592653589793
#define icom (complex<double>(0, 1.))

using namespace std;

//////// INITIALIZATION /////////

double Gb_I, Gb_FS, Gb_SF, R0A,  Ksi_F, Ksi_S, Del0, Del0a, Del0b, Xi1, Xi2, T, w, w1, E1, E, ro_F, ro_S,H, alphaT, aGIP, aGIPlast;
double L_F, L_SL, L_SR, L_S, L_s, L_F1;
int  w_obrez, iter, MODE;
int N_SL, N_SR, N_Mid, N, N_S, N_CPR, NUM_Tech, intG, idGm, idGm1;
double epsG, epsDel, alpha;
double x;
double *Hi, *Ksii, *Roi, *Li, *Stype, *Tci, *Rbi, *a1coef, k1, k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14;
complex<double> *difG, *difGm;

int main()
{   

///////// VARIABLES //////////////

	// Number of Layers

	NUM_Tech=2;
	N=1003;//int(1000*(3+0.015))+20; //cout<<N;
    //   Number of points in grid for each layer    
    //N_Mid=30;// in F layer
    //N_S=
    //N=N_Mid*NUM_Tech;
	N_CPR=1;
    Hi = new double[NUM_Tech];
    Ksii = new double[NUM_Tech]; // Normalized on S
    Roi = new double[NUM_Tech];  // Normalized on S
	Li = new double[NUM_Tech]; 
	Stype = new double[NUM_Tech]; 
	Tci = new double[NUM_Tech];
	Rbi = new double[NUM_Tech];
	difG=new complex<double>[N];
	difGm=new complex<double>[N];
	a1coef=new double[N];
	// Initialization of arrays
    complex<double> *G, *GG, *G1, *GG1, *F,*FF,*Gm, *GGm, *G1m, *GG1m, *Fm,*FFm, *Del, *dDel, *reDelOUT, *imDelOUT;
    double *Is, *Is1;
    G = new complex<double>[N];
    GG = new complex<double>[N];
    GG1 = new complex<double>[N];
    G1 = new complex<double>[N];
	F = new complex<double>[N];
	FF = new complex<double>[N];
	Gm = new complex<double>[N];
    GGm = new complex<double>[N];
    GG1m = new complex<double>[N];
    G1m = new complex<double>[N];
	Fm = new complex<double>[N];
	FFm = new complex<double>[N];
    Del = new complex<double>[N];
    reDelOUT = new complex<double>[N];
    imDelOUT = new complex<double>[N];
    Is = new double[N];
    Is1 = new double[N];
    double** DOSbuf = new double *[4];
		   for (int i = 0; i < 4; i++)
			   DOSbuf[i] = new double[400000];
	double** DOSbuf2 = new double *[4];
		   for (int i = 0; i < 4; i++)
			   DOSbuf2[i] = new double[400000];
	for (int i = 0; i < N; i++)
	{
		G[i] = 1.+0.1*icom;
		G1[i] = 1.+0.1*icom;
		GG1[i] = 1.+0.1*icom;
		GG[i] = 1.+0.1*icom;
		Gm[i] = 1.+0.1*icom;
		G1m[i] = 1.+0.1*icom;
		GG1m[i] = 1.+0.1*icom;
		GGm[i] = 1.+0.1*icom;
	}
//   Size of the layers    
L_S=3; // they are can be changed below in the loops
L_F=0.15;  T=0.05;
    //   Iteration constants    
    epsG=1e-6;  // accuracy of normal Green function G iteration loop
    epsDel=1e-3;  // accuracy of pair potential \Delta iteration loop
    alpha=0.7; // parameter in Delta loop. Just for increase of convergence
    iter=0; // variable to count number of iterations in Delta loop
    MODE=0;
    //   Material parameters      
 Ksi_S=1; // Coherence length of superconductor. Used for normalization of all other lengths and scales
 Ksi_F=1.48; // Coherence length of ferromagnet
 Xi2=0; // Phases of pair potential on the left (Xi1) and right (Xi2) electrodes 
 Xi1=-Xi2; // Josephson phase is their difference: Xi = Xi2 - Xi1
 w_obrez=int(50/T); // Number of Matsubara frequencies used in calculation. 
 H =5.; // Exchange field
 ro_S=520;  // Resistivity of superconductor
 ro_F=250;  // Resistivity of ferromagnet
Gb_FS=3; // boundary parameter of Ferromagnet/Superconductor interface. It depends from resistivity of the interface.
R0A = Gb_FS *ro_F* Ksi_F;
Gb_SF= Gb_FS *ro_F* Ksi_F/ ro_S/ Ksi_S; // it also deterimines Superconductor/Ferromagnet interface
cout<<sqrt(3.-4.*icom)<<endl;
//_____________  Creation of OUTPUT file_______________________//
		string name1("H="),name2("delH="),name3("11for2DF_FH="),str1,str2,str3,str4,str5,str6;// cheate one big file
		stringstream s1,s2,s3,s4,s5,s6;
		s1<<H; s2<<T; s3<<L_S; s4<<L_F; s5<<L_s; s6<<Gb_FS;
		s1>>str1; s2>>str2; s3>>str3; s4>>str4; s5>>str5; s6>>str6;
		name1.append(str1); name1.append("T=");name1.append(str2);name1.append("LS=");name1.append(str3);name1.append("LF=");name1.append(str4);name1.append("gb=");name1.append(str6);name1.append(".txt");
		name2.append(str1); name2.append("T=");name2.append(str2);name2.append("LS=");name2.append(str3);name2.append("LF=");name2.append(str4);name2.append("gb=");name2.append(str6);name2.append(".txt");
		name3.append(str1); name3.append("T=");name3.append(str2);name3.append("LS=");name3.append(str3);name3.append("LF=");name3.append(str4);name3.append("gb=");name3.append(str6);name3.append(".txt");
		const char *file1=name1.c_str();                                                                                                                                                       
		const char *file2=name2.c_str();
		const char *file3=name3.c_str();
		ofstream fout1(file1);
		ofstream fout2(file2);
		ofstream fout3(file3);
	R0A = Gb_FS *ro_F* Ksi_F;
	Gb_SF= Gb_FS *ro_F* Ksi_F/ ro_S/ Ksi_S;
		
/*for(L_F=0.0001;L_F<100001;L_F+=0.001)
if((abs(L_F-0.003)<1e-5))//||(abs(L_F-0.6)<1e-5)||(abs(L_F-1)<1e-5)||(abs(L_F-2)<1e-5))	
{
	//fout1<<"L_F="<<L_F<<endl;
for(Gb_FS=0.01;Gb_FS<100001;Gb_FS+=0.001)
if((abs(Gb_FS-3)<1e-5))//||(abs(Gb_FS-10)<1e-5)||(abs(Gb_FS-30)<1e-5))//||(abs(Gb_FS-1)<1e-5)||(abs(Gb_FS-5)<1e-5)||(abs(Gb_FS-10)<1e-5))
//if((abs(Gb_FS-0.25)<1e-5)||(abs(H-5)<1e-5)||(abs(H-10)<1e-5)||(abs(H-15)<1e-5)||(abs(H-20)<1e-5))
{
	//fout1<<"Gb_FS="<<Gb_FS<<endl;/*
for(alphaT=0;alphaT<10.1;alphaT+=0.001)
if((abs(alphaT-0.25)<1e-5)||(abs(alphaT-0.25)<1e-5)||(abs(alphaT-0.5)<1e-5)||(abs(alphaT-1)<1e-5))//||(abs(alphaT-2)<1e-5)||(abs(alphaT-5)<1e-5)||(abs(alphaT-10)<1e-5))
{/*/

//	fout1<<"alphaT="<<alphaT<<endl;
	
	for (int i = 0; i < NUM_Tech; i++)
	{
		if (i == 0)
		{
			 Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S; Stype[i] = 1.; Tci[i] =1; Rbi[i] = R0A;
		}
		if (i == 1)
		{
			Hi[i] = H; Ksii[i] = Ksi_F; Roi[i] = ro_F; Li[i] = L_F; Stype[i] = 0.;  Tci[i] = 0.; Rbi[i] = R0A;
		}
	}           
 	N_Mid=101;//int(200*L_F)+1;;// in F layer
    N_S=1;//int(L_S/get_h(N_S))+1;
    N=N_Mid+N_S; cout<<N_Mid<<"	"<<N_S<<"	"<<N<<"	"<<get_h(N_S-1)<<"	"<<get_h(N_S)<<endl;
Del0 = SelfConsZero();// Calculation of initial Delta for iterations. It is pair potential of bulk material for certain temperature.
 	complex<double> x1,x2,x3,x4,x5,x6;
 	ifstream file ("C:\\Новая папка (2)\\учёбка\\научная деятельность\\задание\\Tc versus df\\эксперимент с анизотропией\\forGOODcalcdelH=5T=0.05LS=3LF=0.15gb=3.txt");
		if (file.is_open()) cout << "good\n";
		else cout << "beed\n\n" << endl;	
		for (int i=0; i<N; i++)      file>>x1>>x2>>reDelOUT[i]>>x3;	for(int i=0; i<N;++i) Del[i]=reDelOUT[i]+0.*icom;	
//SelfCons(G, Del, 0);     // Self-consistent  calculation of pair potential.
	for (int i=0; i<N_S; i++)      	Del[i]=Del0;    for (int i=N_S; i<N; i++)   Del[i]=0;
   ofstream fout("G0 for H=5 del0.txt"); int i=0; 
    Gcalc0(N_S,w,Del);
	//for (Hi[1]=0.0; Hi[1]<10; Hi[1]+=0.01)
	//if(abs(Hi[1]-5)<1e-5)//||(abs(Hi[1]-pi/16)<1e-5)||(abs(Hi[1]-pi/8)<1e-5)||(abs(Hi[1]-pi/4)<1e-5)||(abs(Hi[1]-pi/2)<1e-5)||(abs(Hi[1]-pi)<1e-5)||(abs(Hi[1]-2*pi)<1e-5))
	/*for (E=-10; E<10; E+=0.01)
	{	complex<double> w = -icom*E ;
		G[i]=Gcalc0(i,w);		x1=Gcalc0(i,w);	
		G1[i]=Gcalc0(i,-w);		x2=Gcalc0(i,-w);  
		Hi[1]*=-1;
		fout<<fixed<<E<<" "<<Hi[1]<<" "<<real(G[i])<<" "<<real(G1[i])<<" "<<imag(G1[i])<<" "<<imag(G[i])<<; 
		G[i]=Gcalc0(i,w);		x3=Gcalc0(i,w);	
		G1[i]=Gcalc0(i,-w); 	x4=Gcalc0(i,-w);
		Hi[1]*=-1;
		fout<<fixed<<" "<<real(G[i])<<" "<<imag(G[i])<<" "<<real(G1[i])<<" "<<imag(G1[i])<<endl;
		if(abs(abs(imag(x1))-abs(imag(x4)))>1e-7) cout<<"bed +w"<<w<<endl;//
		if(abs(abs(imag(x2))-abs(imag(x3)))>1e-7) cout<<"bed -w"<<w<<endl;
	}//*/
	i=N-1;
	//for(Hi[1]=5; Hi[1]>2; Hi[1]-=0.5)
	
	
	for(k1=-1.;k1<1.9; k1+=2.)
	for(k2=-1.;k2<1.9; k2+=2.)
	for(k3=-1.;k3<1.9; k3+=2.)
	for(k4=-1.;k4<1.9; k4+=2.)
	for(k5=-1.;k5<1.9; k5+=2.)
	for(k6=-1.;k6<1.9; k6+=2.)
	for(k7=-1.;k7<1.9; k7+=1.)
	for(k8=-1.;k8<1.9; k8+=2.)
	for(k9=-1.;k9<1.9; k9+=2.)
	for(k10=-1.;k10<1.9; k10+=1.)
	for(k11=-1.;k11<1.9; k11+=1.)
	for(k12=-1.;k12<1.9; k12+=2.)
	{	i=1;
		for (E=-10; E<10; E+=0.1)
		{	complex<double> w = -icom*E ;
			x1=Gcalc0(0,w,Del);	x2=Gcalc0(0,-w,Del); x3=Gcalc0(N_S,w,Del); x4=Gcalc0(N_S,-w,Del); Hi[1]*=-1; x5=Gcalc0(N_S,w,Del); x6=Gcalc0(N_S,-w,Del); Hi[1]*=-1;
			//fout<<fixed<<E<<" "<<k1<<" "<<k2<<" "<<k3<<" "<<k4<<" "<<k5<<" "<<k6<<" "<<k7<<"  "<<Hi[1]<<" "<<real(x1)<<" "<<real(x2)<<" "<<imag(x1)<<" "<<imag(x2)<<" "<<" "<<real(x3)<<" "<<real(x4)<<" "<<real(x5)<<" "<<real(x6)<<" "<<imag(x3)<<" "<<imag(x4)<<" "<<imag(x5)<<" "<<imag(x6)<<endl;//" "<<SQRT(0, w, Del,0.)<<" "<<SQRT(0, -w, Del,0.)<<endl;
			/*if((abs(abs(imag(x1))-abs(imag(x2)))>1e-7)||(abs((real(x1))-(real(x2)))>1e-9)||(abs(imag(x1)+imag(x2))>1e-7)) cout<<"bed +w in S"<<w<<"	"<<abs(imag(x1)+imag(x2))<<endl;//
			if((abs(abs(imag(x3))-abs(imag(x6)))>1e-7)||(abs((real(x3))-(real(x6)))>1e-9)||(abs(imag(x3)+imag(x6))>1e-7)) cout<<"bed +w in F"<<w<<"	"<<abs(imag(x3)+imag(x6))<<endl;//
			if((abs(abs(imag(x4))-abs(imag(x5)))>1e-7)||(abs((real(x4))-(real(x5)))>1e-9)||(abs(imag(x4)+imag(x5))>1e-7)) cout<<"bed -w in F"<<w<<"	"<<abs(imag(x4)+imag(x5))<<endl;
			/*/
			if((real(x3)<-1e-10)||(real(x4)<-1e-10)||(real(x5)<-1e-10)||(real(x6)<-1e-10)||(abs(abs(imag(x3))-abs(imag(x6)))>1e-7)||(abs((real(x3))-(real(x6)))>1e-9)||(abs(imag(x3)+imag(x6))>1e-7)||(abs(abs(imag(x4))-abs(imag(x5)))>1e-7)||(abs((real(x4))-(real(x5)))>1e-9)||(abs(imag(x4)+imag(x5))>1e-7)) 
			i=0;
		}
		if(i==1) cout<<k1<<" "<<k2<<" "<<k3<<" "<<k4<<" "<<k5<<" "<<k6<<" "<<k7<<" "<<k8<<" "<<k9<<" "<<k10<<" "<<k11<<endl;
	}
	//*/
//____________main calculate________________________________________________
	for(int n=0; n<N;++n) //remember Del(x}
	fout2<<fixed<<n<<"\t" <<H <<"\t" <<  real(Del[n])<<"\t" << imag(Del[n])<<endl;		
    	double dE=0.025, IG1=0,IG2=0,IG3=0,G11,G2,G3,G4; int n=0; 
		aGIPlast=  0.02500;//aGIPlast=aGIP; 
		cout<<"w-delEimag(w)/abs(w"<<endl;
		
		//for(aGIP=1.500;aGIP>1.400;aGIP-=0.02)
		for (E=-0.4; E<-0.439; E+=dE)
		//if((abs(E+0.4)<1e-5))//||(abs(E+10)<1e-5)||((abs(E)<7)&&(abs(E)>1.749))||((abs(E)<0.48)&&(abs(E)>0.35))||((abs(E)<0.2)&&(abs(E)>0.02)))//&&(abs(E)>1.749))(abs(E+0.4)<1e-5)||(abs(E+0.1)<1e-5))
		{	//double dGG=1; aGIP=0.15;
		//	while(((dGG>1e-4)||(dGG<1e-10))&&(aGIP>0))			{
				//cout<<aGIP<<endl;
				DOS(E,G,G1,F,Fm,Del);
				Hi[1]*=-1;
				//DOS(E,GG,GG1,FF,FFm,Del);
				Hi[1]*=-1;//*/
				//dGG=abs(real(G[N-1])-real(GG[N-1]));
				
		//	}
			//fout1 <<fixed<<E  <<"\t"<<aGIP<<"\t"<<  real(G[10]+GG[10])/2. <<"\t"<< imag(G[10]+GG[10])/2.<<"\t"<<  real(G1[10]+GG1[10])/2. <<"\t"<< imag(G1[10]+GG1[10])/2.<<"\t"<<  real(G[N-1]+GG[N-1])/2. <<"\t"<< imag(G[N-1]+GG[N-1])/2.<<"\t"<<  real(G1[N-1]+GG1[N-1])/2. <<"\t"<< imag(G1[N-1]+GG1[N-1])/2.<<endl;//*/
			fout1 <<fixed<<E  <<"  " <<aGIP <<"  "<< real(G[10]) <<" "<< imag(G[10])<<" "<< real(G1[10]) <<" "<< imag(G1[10])<<"  "<< real(GG[10]) <<" "<< imag(GG[10])<<" "<< real(GG1[10]) <<" "<< imag(GG1[10])<<"  "<< real(G[N-1]) <<" "<< imag(G[N-1])<<" "<< real(G1[N-1]) <<" "<< imag(G1[N-1])<<"  "<< real(GG[N-1]) <<" "<< imag(GG[N-1])<<" "<< real(GG1[N-1]) <<" "<< imag(GG1[N-1])<<endl;
			 
				//for (int i=N-1; i<N; i++)		fout3 <<fixed<<E << real(G[i]) <<" "<< imag(G[i])<<" "<< real(G1[i]) <<" "<< imag(G1[i])<<" "<< real(GG[i]) <<" "<< imag(GG[i])<<" "<< real(GG1[i]) <<" "<< imag(GG1[i])<<endl;
			 /*for(aGIP=0.1; aGIP<1.2; aGIP+=0.1)
			 {DOS(E,G,G1,F,Fm,Del); fout1 <<fixed<<E  <<"  " <<aGIP <<"  "<< real(G[110]) <<" "<< imag(G[110])<<" "<< real(G1[110]) <<" "<< imag(G1[110])<<"  "<< real(G[N-1]) <<" "<< imag(G[N-1])<<" "<< real(G1[N-1]) <<" "<< imag(G1[N-1])<<endl;}
			 for(aGIP=1.2; aGIP<10; aGIP+=0.5)
			 {DOS(E,G,G1,F,Fm,Del); fout1 <<fixed<<E  <<"  " <<aGIP <<"  "<< real(G[110]) <<" "<< imag(G[110])<<" "<< real(G1[110]) <<" "<< imag(G1[110])<<"  "<< real(G[N-1]) <<" "<< imag(G[N-1])<<" "<< real(G1[N-1]) <<" "<< imag(G1[N-1])<<endl;}//*
			*
			for(aGIP=-10; aGIP<-1.2; aGIP+=1)
			{DOS(E,G,G1,F,Fm,Del);fout1 <<fixed<<E  <<"  " <<aGIP <<"  "<< real(G[110]) <<" "<< imag(G[110])<<" "<< real(G1[110]) <<" "<< imag(G1[110])<<"  "<< real(G[N-1]) <<" "<< imag(G[N-1])<<" "<< real(G1[N-1]) <<" "<< imag(G1[N-1])<<endl;}//*	
			/*/
			//for(aGIP=-0.4; aGIP<0; aGIP+=0.05)//*
			//{DOS(E,G,G1,F,Fm,Del); fout1 <<fixed<<E  <<"  " <<aGIP <<"  "<< real(G[110]) <<" "<< imag(G[110])<<" "<< real(G1[110]) <<" "<< imag(G1[110])<<"  "<< real(G[N-1]) <<" "<< imag(G[N-1])<<" "<< real(G1[N-1]) <<" "<< imag(G1[N-1])<<endl;}	
			
			
			
			for (int i=0; i<N; i++)			fout3 <<fixed<<E  <<" "<<i<<" "<< real(F[i]) <<" "<< imag(F[i])<<" "<< real(Fm[i]) <<" "<< imag(Fm[i])<<" "<< real(FF[i]) <<" "<< imag(FF[i])<<" "<< real(FFm[i]) <<" "<< imag(FFm[i])<<" "<< real(G[i]) <<" "<< imag(G[i])<<" "<< real(G1[i]) <<" "<< imag(G1[i])<<" "<< real(GG[i]) <<" "<< imag(GG[i])<<" "<< real(GG1[i]) <<" "<< imag(GG1[i])<<endl;
			//fout3 <<fixed<<E  <<" "<<i<<" "<< real(F[i]) <<" "<< imag(F[i])<<" "<< real(Fm[i]) <<" "<< imag(Fm[i])<<" "<< real(G[i]) <<" "<< imag(G[i])<<" "<< real(G1[i]) <<" "<< imag(G1[i])<<endl;
			//fout1 <<fixed<<E  <<"  "<< real(G[110]) <<" "<< imag(G[110])<<" "<< real(G1[110]) <<" "<< imag(G1[110])<<"  "<< real(G[N-1]) <<" "<< imag(G[N-1])<<" "<< real(G1[N-1]) <<" "<< imag(G1[N-1])<<endl;
			// fout3 <<fixed<<E  <<"\t"<<aGIP<<"\t"<<i<<"\t"<< F[i] <<"\t"<< Fm[i]<<"\t" << G[i] <<"\t"<< G1[i]<<"\t"<< FF[i] <<"\t"<< FFm[i]<<"\t" << GG[i] <<"\t"<< GG1[i]<<endl;	
			//for (int i=0; i<N; i++) fout3 <<fixed<<E  <<"\t"<<i<<"\t"<< real(F[i]+FF[i])/2. <<"\t"<< imag(F[i]+FF[i])/2.<<"\t"<< real(Fm[i]+FFm[i])/2. <<"\t"<< imag(Fm[i]+FFm[i])/2.<<"\t" << real(G[i]+GG[i])/2. <<"\t"<< imag(G[i]+GG[i])/2. <<"\t" << real(G1[i]+GG1[i])/2. <<"\t"<< imag(G1[i]+GG1[i])/2.<<endl;
			
			//fout1 <<fixed<< E << "\t"<<alphaT<< "\t" <<  real(G[N-1])<< endl;
			//fout1 <<fixed<< E << "\t" <<  2.*real(G[N-1])/2.<< endl;//for (int i=0; i<N; i++)	fout1 <<fixed<<E  << "\t"<< i<<"\t"<< real(F[i]*G[i]/(-icom*E+1e-5)*conj(F[i]*G[i]/(-icom*E+1e-5)))<<"\t"<<imag(F[i]*G[i]/(-icom*E+1e-5)*conj(F[i]*G[i]/(-icom*E+1e-5)))<<"\t"<< real(G[i]*G[i])<<"\t"<< imag(G[i]*G[i])<<"\t"<<(-F[i]*G[i]/(-icom*E+1e-5)*(F[i]*G[i]/(-icom*E+1e-5))+G[i]*G[i])<<endl;
			
			/*DOSbuf[0][n]=real(G[0]); DOSbuf[1][n]=real(G[99]); DOSbuf[2][n]=real(G[100]); DOSbuf[3][n]=real(G[199]);
			if(abs(E)<1e-9){ DOSbuf[0][n]=DOSbuf[0][n-1]; DOSbuf[1][n]=DOSbuf[1][n-1];DOSbuf[2][n]=DOSbuf[2][n-1];DOSbuf[3][n]=DOSbuf[3][n-1];}	real(F[i]*G[i]/sqrt(icom*E*icom*E)+F[i]*conj(Fm[i]))
			n++;
			/*if(abs(E+dE)<1e-9){ G11=real(G[1]); G2=real(G[99]); G3=real(G[100]); G4=real(G[199]);}
			if(abs(E-0)>1e-9)
				fout1 <<fixed<< E << "\t"<<H<< "\t"<< real(G[0]) << "\t" << real(G[99]) <<"\t" << real(G[100])<<"\t" << real(G[199])<< endl;
			else 
				fout1 <<fixed<< E << "\t"<<H<< "\t"<< G11 << "\t" << G2 <<"\t" << G3 <<"\t" << G4 << "\t"<<endl;//*/
		}
		/*for (E=6; E<200; E+=dE)
		{
			DOS(E,G,G1,F,Del);
			DOSbuf[0][n]=real(G[0]); DOSbuf[1][n]=real(G[99]); DOSbuf[2][n]=real(G[100]); DOSbuf[3][n]=real(G[199]);
			if(abs(E)<1e-9){ DOSbuf[0][n]=DOSbuf[0][n-1]; DOSbuf[1][n]=DOSbuf[1][n-1];DOSbuf[2][n]=DOSbuf[2][n-1];DOSbuf[3][n]=DOSbuf[3][n-1];}
			n++;
			/*if(abs(E+dE)<1e-9){ G11=real(G[1]); G2=real(G[99]); G3=real(G[100]); G4=real(G[199]);}
			if(abs(E-0)>1e-9)
				fout1 <<fixed<< E << "\t"<<H<< "\t"<< real(G[0]) << "\t" << real(G[99]) <<"\t" << real(G[100])<<"\t" << real(G[199])<< endl;
			else 
				fout1 <<fixed<< E << "\t"<<H<< "\t"<< G11 << "\t" << G2 <<"\t" << G3 <<"\t" << G4 << "\t"<<endl;//*
		}//*
		
	Hi[1]*=-1;
		n=0; 
		for (E=-6; E<6; E+=dE)
		{
			DOS(E,G,G1,F,Del);
			DOSbuf2[0][n]=real(G[0]); DOSbuf2[1][n]=real(G[99]); DOSbuf2[2][n]=real(G[100]); DOSbuf2[3][n]=real(G[199]);
			if(abs(E)<1e-9) { DOSbuf2[0][n]=DOSbuf2[0][n-1]; DOSbuf2[1][n]=DOSbuf2[1][n-1];DOSbuf2[2][n]=DOSbuf2[2][n-1];DOSbuf2[3][n]=DOSbuf2[3][n-1];}
			n++;
			/*if(abs(E+dE)<1e-9){ G11=real(G[1]); G2=real(G[99]); G3=real(G[100]); G4=real(G[199]);}
			if(abs(E-0)>1e-9)
				//fout1 <<fixed<< E << "\t"<<H<< "\t"<< real(G[0]) << "\t" << real(G[99]) <<"\t" << real(G[100])<<"\t" << real(G[199])<< endl;
			else 
				//fout1 <<fixed<< E << "\t"<<H<< "\t"<< G11 << "\t" << G2 <<"\t" << G3 <<"\t" << G4 << "\t"<<endl;/*
		}
		/*
		for (E=6; E<200; E+=dE)
		{
			DOS(E,G,G1,F,Del);
			DOSbuf2[0][n]=real(G[0]); DOSbuf2[1][n]=real(G[99]); DOSbuf2[2][n]=real(G[100]); DOSbuf2[3][n]=real(G[199]);
			if(abs(E)<1e-9) { DOSbuf2[0][n]=DOSbuf2[0][n-1]; DOSbuf2[1][n]=DOSbuf2[1][n-1];DOSbuf2[2][n]=DOSbuf2[2][n-1];DOSbuf2[3][n]=DOSbuf2[3][n-1];}
			n++;
			/*if(abs(E+dE)<1e-9){ G11=real(G[1]); G2=real(G[99]); G3=real(G[100]); G4=real(G[199]);}
			if(abs(E-0)>1e-9)
				//fout1 <<fixed<< E << "\t"<<H<< "\t"<< real(G[0]) << "\t" << real(G[99]) <<"\t" << real(G[100])<<"\t" << real(G[199])<< endl;
			else 
				//fout1 <<fixed<< E << "\t"<<H<< "\t"<< G11 << "\t" << G2 <<"\t" << G3 <<"\t" << G4 << "\t"<<endl;//*
		}	//*/
	
	/*for(int ii=0;ii<(200-6)*100;ii++)
	{	
			fout1 <<fixed<< ((ii-20000)*0.01) << "\t" <<abs(H)<< "\t"<< ((DOSbuf[0][ii]+DOSbuf2[0][ii])/2.)<< "\t" <<((DOSbuf[1][ii]+DOSbuf2[1][ii])/2.)<< "\t" <<((DOSbuf[2][ii]+DOSbuf2[2][ii])/2.)<< "\t" <<((DOSbuf[3][ii]+DOSbuf2[3][ii])/2.)<< "\t" <<endl;
			
	}
	for(int ii=(200-6)*100;ii<(400-12)*100;ii++)
	{	
			fout1 <<fixed<< ((ii-20000+1200)*0.01) << "\t" <<abs(H)<< "\t"<< ((DOSbuf[0][ii]+DOSbuf2[0][ii])/2.)<< "\t" <<((DOSbuf[1][ii]+DOSbuf2[1][ii])/2.)<< "\t" <<((DOSbuf[2][ii]+DOSbuf2[2][ii])/2.)<< "\t" <<((DOSbuf[3][ii]+DOSbuf2[3][ii])/2.)<< "\t" <<endl;
			
	}//*

	for(int ii=(20000-600);ii<(40000-1200);ii++)
	{
			fout1 <<fixed<< (ii+1200-20000) << "\t" <<abs(H)<< "\t"<< ((DOSbuf[0][ii]+DOSbuf2[0][ii])/2.)<< "\t" <<((DOSbuf[1][ii]+DOSbuf2[1][ii])/2.)<< "\t" <<((DOSbuf[2][ii]+DOSbuf2[2][ii])/2.)<< "\t" <<((DOSbuf[3][ii]+DOSbuf2[3][ii])/2.)<< "\t" <<endl;
	}//*/
//}}}
		
		//fout1<<"\n"<<"H="<<H<<", T="<<T<<", L_S="<<L_S<<", L_F="<<L_F<<", Gb_FS="<<Gb_FS<<", ro_S="<<ro_S<<", ro_F="<<ro_F<<", R0A="<<R0A<<"\n";
		//fout2<<"\n"<<"H="<<H<<", T="<<T<<", L_S="<<L_S<<", L_F="<<L_F<<", Gb_FS="<<Gb_FS<<", ro_S="<<ro_S<<", ro_F="<<ro_F<<", R0A="<<R0A<<"\n";
   
	
    system("PAUSE");
    return EXIT_SUCCESS;
}
