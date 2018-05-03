#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>

void charge(Int_t run){

//=====Run number and scaler of fast bus

   //Int_t run, I;
   //cout << "\n Please enter a Run Number (-1 to exit):  ";
   //cin >>run;
   //if(run==-1) return;

   //cout << "\n Please enter the asked current (-1 to exit) ";
   //cin >>I;
   //if(I==-1) return;

//============  Reading the Rootfile =======================//

  const TString rootfilePath = "/Users/nathalysantiesteban/Documents/JLabWork/Tritium/";
  std::ostringstream str;
  str << rootfilePath<<"tritium_"<<run;
  TString basename = str.str().c_str();
  TString rootfile = basename + ".root";
  TChain* T;
  TChain*ev;

  if(run < 20000){
    T = new TChain("TSLeft");
    ev = new TChain("evLeft");
  }
  if(run > 20000){
    T = new TChain("TSRight");
    ev = new TChain("evRight");
  }

  Long_t split=0;

  char* file = 0;

  //====adding splits rootfiles =======================//

  Long_t i=0;
  while ( !gSystem->AccessPathName(rootfile.Data()) ) {
     T->Add(rootfile.Data());
     ev->Add(rootfile.Data());
     cout << "ROOT file " << rootfile << " added to TChain." << endl;
     i++;
     rootfile = basename + "_" + i + ".root";
   }

  if(!T->GetEntries()){
     cerr<< "Not root file was found" << endl;
     return;
  }

    //==finish adding splits rootfiles=====================//  

  Double_t d3c, d10c, unewc, dnewc;
  Double_t d3r, d10r, unewr, dnewr;
  Double_t clk;

  ev->SetBranchAddress("evLeftd3",&d3c);
  ev->SetBranchAddress("evLeftd10",&d10c);
  ev->SetBranchAddress("evLeftunew",&unewc);
  ev->SetBranchAddress("evLeftdnew",&dnewc);

  ev->SetBranchAddress("evLeftd3_r",&d3r);
  ev->SetBranchAddress("evLeftd10_r",&d10r);
  ev->SetBranchAddress("evLeftunew_r",&unewr);
  ev->SetBranchAddress("evLeftdnew_r",&dnewr);

  ev->SetBranchAddress("evLeftLclock",&clk);
  
  //========Calibration Values of the bcm accorfing to the bcm[] array

   Double_t g[4] = {0.0001071,3.795e-5,0.0003975,0.0003281};
   Double_t ger[4] = {4.481e-7,4.087e-7,1.663e-6,1.32e-6};

   //Double_t of[4] = {0,0,0,0};
   Double_t of[4] = {0.2659,0.007335,0.1077,0.1402};
   Double_t ofer[4] = {0.1027, 0.1464, 0.1402, 0.1033};

    //============== Variables ====================================================//


   Int_t evnentries = ev->GetEntries();

   //==========preliminary charge with the final and initial values ======//
  
  Double_t cd3_i, cd10_i, cunew_i, cdnew_i, clk_i, cd3_f, cd10_f, cunew_f, cdnew_f, clk_f;
  Double_t Qd3_tot, Qd10_tot, Qunew_tot, Qdnew_tot;

  ev->GetEntry(0);
  cd3_i = d3c; cd10_i = d10c; cunew_i = unewc; cdnew_i = dnewc;
  clk_i = clk;

  ev->GetEntry(evnentries-1);
  cd3_f = d3c; cd10_f = d10c; cunew_f = unewc; cdnew_f = dnewc;
  clk_f = clk;
  
  Qd3_tot = (g[0]*(cd3_f-cd3_i))+ ((of[0]*(clk_f-clk_i))/103700.0 );
  Qd10_tot = (g[1]*(cd10_f-cd10_i))+ ((of[1]*(clk_f-clk_i))/103700.0 );
  Qunew_tot = (g[2]*(cunew_f-cunew_i))+ ((of[2]*(clk_f-clk_i))/103700.0 );
  Qdnew_tot = (g[3]*(cdnew_f-cdnew_i))+ ((of[3]*(clk_f-clk_i))/103700.0 );

  cout << "The preliminary charge for d10: " << Qd3_tot << " uC" << endl;
  cout << "The preliminary charge for d10: " << Qd10_tot << " uC" << endl;
  cout << "The preliminary charge for unew: " << Qunew_tot << " uC" << endl;
  cout << "The preliminary charge for dnew: " << Qdnew_tot << " uC" << endl;

   //==========Current Calculation==============================================

  Double_t I_d3[evnentries], I_d10[evnentries], I_unew[evnentries], I_dnew[evnentries];
  Double_t I_d3er[evnentries], I_d10er[evnentries], I_unewer[evnentries], I_dnewer[evnentries];
  
  for (Int_t i = 0; i<evnentries-1 ; i++) {
    ev->GetEntry(i);
    //========current
    I_d3[i] = g[0]*d3r + of[0];
    I_d10[i] = g[1]*d10r + of[1];
    I_unew[i] = g[2]*unewr + of[2];
    I_dnew[i] = g[3]*dnewr + of[3];
    //==========current error from fit errors=====//////////
    I_d3er[i] = sqrt( pow(ger[0]*d3r,2) + pow(ofer[0],2));
    I_d10er[i] = sqrt( pow(ger[1]*d10r,2) + pow(ofer[1],2));
    I_unewer[i] = sqrt( pow(ger[2]*unewr,2) + pow(ofer[2],2));
    I_dnewer[i] = sqrt( pow(ger[3]*dnewr,2) + pow(ofer[3],2));    
  }

   //============   Current Average  ================================================

   Double_t curd3, curd10, curdnew, curunew, co3, co10, counew, codnew;
   Double_t curd3er, curd10er, curdnewer, curunewer;
   Double_t curd3er1, curd10er1, curdnewer1, curunewer1;
   curd3 = curd10 = curdnew = curunew = co3 = co10 = counew = codnew = 0;
   curd3er = curd10er = curdnewer = curunewer = 0; 



   for(Int_t i=0; i<evnentries-1; i++){
    if (I_d3[i]>(0.5) && I_d3[i]<(30)){curd3+=I_d3[i]; co3+=1; }
    if (I_d10[i]>(0.5)&& I_d10[i]<(30)){curd10+=I_d10[i]; co10+=1; }
    if (I_unew[i]>(0.5)&& I_unew[i]<(30)){curunew+=I_unew[i]; counew+=1;}
    if (I_dnew[i]>(0.5)&& I_dnew[i]<(30)){curdnew+=I_dnew[i]; codnew+=1; }
   }
  
  
  curd3/=co3; curd10/=co10; curunew/=counew; curdnew/=codnew;

   


  
   for(Int_t i=0; i<evnentries-1; i++){
    if (I_d3[i]>0.5 && I_d3[i]<(30)){ curd3er +=pow(curd3-I_d3[i],2); curd3er1 += pow(I_d3er[i],2); }
    if (I_d10[i]>0.5 && I_d10[i]<(30)){ curd10er +=pow(curd10-I_d10[i],2); curd10er1 += pow(I_d10er[i],2);}
    if (I_unew[i]>0.5 && I_unew[i]<(30)){curunewer +=pow(curunew-I_unew[i],2); curunewer1 += pow(I_unewer[i],2);}
    if (I_dnew[i]>0.5 && I_dnew[i]<(30)){curdnewer +=pow(curdnew-I_dnew[i],2); curdnewer1 += pow(I_dnewer[i],2);}
   }

   curd3er1 /=co3; curd10er1 /=co10; curunewer1 /=counew; curdnewer1 /=codnew; 

   curd3er = sqrt(curd3er)/co3; curd3er = sqrt(pow(curd3er,2)+curd3er1);
   curd10er = sqrt(curd10er)/co10; curd10er = sqrt(pow(curd10er,2)+curd10er1);
   curunewer = sqrt(curunewer)/counew; curunewer = sqrt(pow(curunewer,2)+curunewer1);
   curdnewer = sqrt(curdnewer)/codnew; curdnewer = sqrt(pow(curdnewer,2)+curdnewer1);
   cout << counew << endl;
   cout << codnew << endl;
 
   cout << curd3 << "  " << curd3er << endl;
   cout << curd10 << "  " << curd10er << endl;
   cout << curunew << "  " << curunewer << endl;
   cout << curdnew << "  " << curdnewer << endl;

   //============   Charge Calculation  ================================================
   
   Double_t Time[evnentries], d3c_t[evnentries], d10c_t[evnentries], unc_t[evnentries], dnc_t[evnentries];
   Double_t d3c_T, d10c_T, unc_T, dnc_T, test;
   d3c_T = d10c_T = unc_T = dnc_T = test=0;

   for (Int_t i = 0; i < evnentries; i++) {
    ev->GetEntry(i);
    Time[i] = clk - test;
    d3c_t[i] = d3c - d3c_T; d10c_t[i] = d10c -  d10c_T;
    unc_t[i] = unewc - unc_T; dnc_t[i] = dnewc - dnc_T;    
    d3c_T = d3c; d10c_T = d10c; unc_T = unewc; dnc_T = dnewc;
    test = clk;
  }


   Double_t  Q_d3[evnentries], Q_d10[evnentries], Q_un[evnentries], Q_dn[evnentries];
   Double_t  Q_d3er[evnentries], Q_d10er[evnentries], Q_uner[evnentries], Q_dner[evnentries];
   
  for (Int_t i = 0; i<evnentries ; i++) {
  //=========charge
    Q_d3[i] = (g[0]*d3c_t[i] + (of[0]*Time[i]/103700.0 ));
    Q_d10[i] = (g[1]*d10c_t[i] + (of[1]*Time[i]/103700.0 ));
    Q_un[i] = (g[2]*unc_t[i] + (of[2]*Time[i]/103700.0 ));
    Q_dn[i] = (g[3]*dnc_t[i] + (of[3]*Time[i]/103700.0 ));

    //===========error of the charge using the fit parameters==========////
    Q_d3er[i] = sqrt( pow(ger[0]*d3c_t[i],2) + pow(ofer[0]*Time[i]/103700.0,2 )); 
    Q_d10er[i] = sqrt( pow(ger[1]*d10c_t[i],2) + pow(ofer[1]*Time[i]/103700.0,2));
    Q_uner[i] = sqrt( pow(ger[2]*unc_t[i],2) + pow(ofer[2]*Time[i]/103700.0,2));
    Q_dner[i] = sqrt( pow(ger[3]*dnc_t[i],2) + pow(ofer[3]*Time[i]/103700.0,2));
  }


  Double_t Charge_d3, Charge_d10, Charge_un, Charge_dn;
  Double_t Charge_d3er, Charge_d10er, Charge_uner, Charge_dner;
  Charge_d3 = Charge_d10 = Charge_un = Charge_dn = 0;
  Charge_d3er = Charge_d10er = Charge_uner = Charge_dner = 0;
  Double_t x[evnentries];

  Double_t f[4];
  Double_t df[4];

  for(Int_t j=0; j<4; j++){
    df[j] = (0.5-of[j])/g[j];
  }

  Int_t tester=0;

  for(Int_t i = 0; i<evnentries-1; i++){
    ev->GetEntry(i);
    x[i] = i; 
    if(d3r >df[0]){ Charge_d3 += Q_d3[i]; Charge_d3er += Q_d3er[i];}
    if(d10r>df[1]){ Charge_d10 += Q_d10[i]; Charge_d10er +=  Q_d10er[i]; }
    if(unewr>df[2]){ Charge_un += Q_un[i]; Charge_uner +=  Q_uner[i]; }
    if(dnewr>df[3]){ Charge_dn += + Q_dn[i]; Charge_dner += Q_dner[i];} 
  }
  
  Charge_d3er = sqrt(Charge_d3er); Charge_d10er = sqrt(Charge_d10er);
  Charge_uner = sqrt(Charge_uner); Charge_dner = sqrt(Charge_dner);

  //Charge_d3er = sqrt( pow((ger[0]*(cd3_f-cd3_i)),2)+ pow(((ofer[0]*(clk_f-clk_i))/103700.0 ),2));
  //Charge_d10er = sqrt( pow((ger[1]*(cd10_f-cd10_i)),2)+ pow(((ofer[1]*(clk_f-clk_i))/103700.0 ),2));
  //Charge_uner = sqrt( pow((ger[2]*(cunew_f-cunew_i)),2)+ pow(((ofer[2]*(clk_f-clk_i))/103700.0 ),2));
  //Charge_dner = sqrt( pow((ger[3]*(cdnew_f-cdnew_i)),2)+ pow(((ofer[3]*(clk_f-clk_i))/103700.0 ),2));


  cout << "Total Charge d3:    " << Charge_d3 << " " << "+/-" << " " << Charge_d3er << " " << "C" << "    or  " << Charge_d3*(1.602e-13) << " " << "+/-" << " " << Charge_d3er*(1.602e-13) << " " << "electrons" << endl;
  cout << "Total Charge d10:    " << Charge_d10 << " " << "+/-" << " " << Charge_d10er << " " << "C" << "  or  " << Charge_d10*(1.602e-13) << " " << "+/-" << " " << Charge_d10er*(1.602e-13) << " " << "electrons" << endl;
  cout << "Total Charge unew:  " << Charge_un << " " << "+/-" << " " << Charge_dner << " " << "C" << "    or  " << Charge_un*(1.602e-13) << " " << "+/-" << " " << Charge_dner*(1.602e-13) << " " << "electrons" <<endl;
  cout << "Total Charge dnew:  " << Charge_dn << " " << "+/-" << " " << Charge_uner << " " << "C" << "    or  " <<  Charge_dn*(1.602e-13) << " " << "+/-" << " " << Charge_uner*(1.602e-13) << " " << "electrons" <<endl;
  
   ofstream outfile ; 
  outfile.open("/Users/nathalysantiesteban/Documents/JLabWork/boiling_tritium/data/charge.txt",ios::in);

  if(!outfile.good()){
    //outfile.close();
    outfile.open("/Users/nathalysantiesteban/Documents/JLabWork/boiling_tritium/data/charge.txt");
    outfile << setiosflags(ios::left) << setw(8) << "run:  "  ;
    outfile << setiosflags(ios::left) << setw(20) << "I(d3) (uA)" ;
    outfile << setiosflags(ios::left) << setw(20) << "I(d3) er(uA)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(d3) (prel uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(d3) (uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(d3) err(uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "I(d10) (uA)" ;
    outfile << setiosflags(ios::left) << setw(20) << "I(d10) er(uA)"; 
    outfile << setiosflags(ios::left) << setw(20) << "Q(d10) (prel uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(d10) (uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(d10) err(uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "I(unew) (uA)" ;
    outfile << setiosflags(ios::left) << setw(20) << "I(unew) er(uA)"; 
    outfile << setiosflags(ios::left) << setw(20) << "Q(unew) (prel uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(unew) (uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(unew) err(uC)";
    outfile << setiosflags(ios::left) << setw(20) << "I(dnew) (uA)" ;
    outfile << setiosflags(ios::left) << setw(20) << "I(dnew) er(uA)"; 
    outfile << setiosflags(ios::left) << setw(20) << "Q(dnew) (prel uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(dnew) (uC)" ;
    outfile << setiosflags(ios::left) << setw(20) << "Q(dnew) err(uC)" << endl;

  }

  outfile.close();
  outfile.open("/Users/nathalysantiesteban/Documents/JLabWork/boiling_tritium/data/charge.txt",ios::out|ios::app);

  //outfile.open("LT_Data.txt");
  outfile<< setiosflags(ios::left) << setw(8) << run ;
  outfile<< setiosflags(ios::left) << setw(20) << curd3;
  outfile<< setiosflags(ios::left) << setw(20) << curd3er;
  outfile<< setiosflags(ios::left) << setw(20) << Qd3_tot;
  outfile << setiosflags(ios::left) << setw(20) << Charge_d3;
  outfile << setiosflags(ios::left) << setw(20) << Charge_d3er;
  outfile<< setiosflags(ios::left) << setw(20) << curd10;
  outfile<< setiosflags(ios::left) << setw(20) << curd10er;
  outfile<< setiosflags(ios::left) << setw(20) << Qd10_tot;
  outfile << setiosflags(ios::left) << setw(20) << Charge_d10;
  outfile << setiosflags(ios::left) << setw(20) << Charge_d10er;
  outfile<< setiosflags(ios::left) << setw(20) << curunew;
  outfile<< setiosflags(ios::left) << setw(20) << curunewer;
  outfile<< setiosflags(ios::left) << setw(20) << Qunew_tot;
  outfile << setiosflags(ios::left) << setw(20) << Charge_un;
  outfile << setiosflags(ios::left) << setw(20) << Charge_uner;
  outfile<< setiosflags(ios::left) << setw(20) << curdnew;
  outfile<< setiosflags(ios::left) << setw(20) << curdnewer;
  outfile<< setiosflags(ios::left) << setw(20) << Qdnew_tot;
  outfile << setiosflags(ios::left) << setw(20) << Charge_dn;
  outfile << setiosflags(ios::left) << setw(20) << Charge_dner << endl;




  outfile.close();



  TCanvas *c1 = new TCanvas("c1","currents",2000,800);
  gStyle->SetOptStat(0);
  c1->Divide(4,2);   

  c1->cd(1);
  TGraphErrors *h3 = new TGraphErrors(evnentries,x,I_d3,0,I_d3er);
  h3->SetTitle("d3 Current");
  h3->GetXaxis()->SetTitle("Entry Number");
  h3->GetYaxis()->SetTitle("Current (uA)");
  h3->GetYaxis()->SetRangeUser(0,25);
  h3->SetMarkerColor(4);
  h3->SetMarkerStyle(7);
  h3->Draw("AP");

  c1->cd(2);
  TGraphErrors *h4 = new TGraphErrors(evnentries,x,I_d10,0,I_d10er);
  h4->SetTitle("d10 Current");
  h4->GetXaxis()->SetTitle("Entry Number");
  h4->GetYaxis()->SetTitle("Current (uA)");
  h4->GetYaxis()->SetRangeUser(0,25);
  h4->SetMarkerColor(4);
  h4->SetMarkerStyle(7);
  h4->Draw("AP");

  c1->cd(3);
  TGraphErrors *h5 = new TGraphErrors(evnentries,x,I_unew,0,I_unewer);
  h5->SetTitle("unew Current");
  h5->GetXaxis()->SetTitle("Entry Number");
  h5->GetYaxis()->SetTitle("Current (uA)");
  h5->GetYaxis()->SetRangeUser(0,25);
  h5->SetMarkerColor(4);
  h5->SetMarkerStyle(7);
  h5->Draw("AP");

  c1->cd(4);
  TGraphErrors *h6 = new TGraphErrors(evnentries,x,I_dnew,0,I_dnewer);
  h6->SetTitle("dnew Current");
  h6->GetXaxis()->SetTitle("Entry Number");
  h6->GetYaxis()->SetTitle("Current (uA)");
  h6->GetYaxis()->SetRangeUser(0,25);
  h6->SetMarkerColor(4);
  h6->SetMarkerStyle(7);
  h6->Draw("AP");
       
  c1->cd(5);
  TGraphErrors *ch3 = new TGraphErrors(evnentries,x,Q_d3,0,Q_d3er);
  ch3->SetTitle("d3 Charge");
  ch3->GetXaxis()->SetTitle("Entry Number");
  ch3->GetYaxis()->SetTitle("Current (uC)");
  ch3->SetMarkerColor(4);
  ch3->SetMarkerStyle(7);
  ch3->Draw("AP");

  c1->cd(6);
  TGraphErrors *ch4 = new TGraphErrors(evnentries,x,Q_d10,0,Q_d10er);
  ch4->SetTitle("d10 Charge");
  ch4->GetXaxis()->SetTitle("Entry Number");
  ch4->GetYaxis()->SetTitle("Current (uC)");
  ch4->SetMarkerColor(4);
  ch4->SetMarkerStyle(7);
  ch4->Draw("AP");

  c1->cd(7);
  TGraphErrors *ch5 = new TGraphErrors(evnentries,x,Q_un,0,Q_uner);
  ch5->SetTitle("un Charge");
  ch5->GetXaxis()->SetTitle("Entry Number");
  ch5->GetYaxis()->SetTitle("Current (uC)");
  ch5->SetMarkerColor(4);
  ch5->SetMarkerStyle(7);
  ch5->Draw("AP");

  c1->cd(8);
  TGraphErrors *ch6 = new TGraphErrors(evnentries,x,Q_dn,0,Q_dner);
  ch6->SetTitle("dn Charge");
  ch6->GetXaxis()->SetTitle("Entry Number");
  ch6->GetYaxis()->SetTitle("Current (uC)");
  ch6->SetMarkerColor(4);
  ch6->SetMarkerStyle(7);
  ch6->Draw("AP");

  //c1->SaveAs(Form("/Users/nathalysantiesteban/Documents/JLabWork/analysis_plots/boiling/charge/Current_run%i_%fuA.pdf",run,I));

}
