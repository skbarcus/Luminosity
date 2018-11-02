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

void BCM_Charges(Int_t run){

//============  Reading the Rootfile =======================//

  const TString rootfilePath = "/home/skbarcus/Tritium/Analysis/He3/Rootfiles/";
  std::ostringstream str;
  str << rootfilePath<<"e08014_"<<run;
  TString basename = str.str().c_str();
  TString rootfile = basename + ".root";
  TChain* T;
  TChain*ev;

  if(run < 20000)     //Must use right scalers as left ones had issues.
    {
      T = new TChain("T");
      ev = new TChain("RIGHT");
    }
  if(run > 20000)
    {
      T = new TChain("T");
      ev = new TChain("EVRIGHT");
    }

  Long_t split=0;

  char* file = 0;

  //====adding splits rootfiles =======================//

  Long_t u=0;
  while ( !gSystem->AccessPathName(rootfile.Data()) ) 
    {
     T->Add(rootfile.Data());
     ev->Add(rootfile.Data());
     cout << "ROOT file " << rootfile << " added to TChain." << endl;
     u++;
     rootfile = basename + "_" + u + ".root";
   }

  if(!T->GetEntries())
    {
     cerr<< "Not root file was found" << endl;
     return;
  }
   //==finish adding splits rootfiles=====================//

  T->SetBranchStatus("*",0);
  T->SetBranchStatus("right_bcm_u1c",1);
  T->SetBranchStatus("right_bcm_d1c",1);
  T->SetBranchStatus("right_clkcount",1);
  ev->SetBranchStatus("*",0);
  ev->SetBranchStatus("bcm_u1c",1);
  ev->SetBranchStatus("bcm_d1c",1);
  ev->SetBranchStatus("clkcount",1);
 

  Int_t edge_off = 2;       //Number of scaler readouts after beam trip waiting for beam to stabilize.  
  Double_t u1c,d1c,rclk;
  Double_t u1c_i, u1c_f, d1c_i, d1c_f, rclk_i,rclk_f;
  Double_t Qu1_rough, Qd1_rough;           //Charges.
  Double_t time_rough = 0.;
  Double_t clk_r = 1024.;      //Clock frequency (Hz).
  Double_t I_set = 120.*0.90;       //Target current for the beam to be set at.
  Double_t e = 1.602E-19;      //e- charge (C).

  //These next two sets of calibration values are from BCM Calibration for E08-014, Patricia Solvignon July 2012.
  /*
  if(run<4000)
    {
      Double_t Cv2f_u1 = 2064./1032., Koff_u1 = 433./217., Cv2f_u3 = 6358./3174., Koff_u3 = -1./66., Cv2f_d1 = 1251./1251., Koff_d1 = 44./42., Cv2f_d3 = 7751./3870., Koff_d3 = -245./-6.;
    }

  if(run>4000)
    {
      Double_t Cv2f_u1 = 2075./1038., Koff_u1 = 366./183., Cv2f_u3 = 6391./3181., Koff_u3 = -206./283., Cv2f_d1 = 1274./1274., Koff_d1 = -3./-5., Cv2f_d3 = 7834./3934., Koff_d3 = -823./49.;
      }*/

  //These calibration values are from Dien which she says came from Zhihong. 
  //Double_t Cv2f_u1 = 2070., Koff_u1 = 400., Cv2f_u3 = 0., Koff_u3 = 0., Cv2f_d1 = 1263., Koff_d1 = 20., Cv2f_d3 = 0., Koff_d3 = 0.;
  //Ok these values actually seem to work.
  Double_t Cv2f_u1 = 1035., Koff_u1 = 200., Cv2f_u3 = 0., Koff_u3 = 0., Cv2f_d1 = 1263., Koff_d1 = 37., Cv2f_d3 = 0., Koff_d3 = 0.;

  /*
  T->SetBranchAddress("right_bcm_u1c",&u1c);
  T->SetBranchAddress("right_bcm_d1c",&d1c);
  T->SetBranchAddress("right_clkcount",&rclk);
  Int_t temp_nevts = T->GetEntries();
  */

  ev->SetBranchAddress("bcm_u1c",&u1c);
  ev->SetBranchAddress("bcm_d1c",&d1c);
  ev->SetBranchAddress("clkcount",&rclk);
  Int_t temp_nevts = ev->GetEntries();
 
  std::ofstream output1 (Form("%d_u1.txt",run), std::ofstream::out);
  std::ofstream output2 (Form("%d_d1.txt",run), std::ofstream::out);

  const int nevts = temp_nevts;    //This recasting of the variable allows me to use it in array lengths!

  cout<<"Number of entries = "<<nevts<<endl;
  
  ev->GetEntry(0);
  u1c_i = u1c;
  d1c_i = d1c;
  rclk_i = rclk;
  cout<<"u1c_i = "<<u1c_i<<"   d1c_i = "<<d1c_i<<"   rclk_i = "<<rclk_i<<endl;

  ev->GetEntry(nevts-1);
  u1c_f = u1c;
  d1c_f = d1c;
  rclk_f = rclk;
  cout<<"u1c_f = "<<u1c_f<<"   d1c_f = "<<d1c_f<<"   rclk_f = "<<rclk_f<<endl;

  time_rough = rclk_f/clk_r;
  cout<<"time_rough = "<<time_rough<<" seconds.   "<<time_rough/60.<<" minutes."<<endl;

  //Calculate a rough beam charge not accounting for beam trips.
  Qu1_rough = ( (u1c_f/time_rough-Koff_u1)/(Cv2f_u1) )*time_rough;
  Qd1_rough = ( (d1c_f/time_rough-Koff_d1)/(Cv2f_d1) )*time_rough;
  cout<<"Qu1_rough = "<<Qu1_rough<<" uC"<<"   Qd1_rough = "<<Qd1_rough<<" uC"<<endl;

  //const int y = nevts;

  Double_t Iu1[nevts], Id1[nevts], u1c_a[nevts], d1c_a[nevts], rclk_a[nevts];
 //Double_t I, u1c_a, rclk_a;


  //Read in all scaler values for the clock and bcms.
  for (Int_t i=0; i<nevts; i++) 
    {
      ev->GetEntry(i);
      u1c_a[i] = u1c;
      d1c_a[i] = d1c;
      rclk_a[i] = rclk;
      //cout<<"u1c_a["<<i<<"] = "<<u1c_a[i]<<"   rclk["<<i<<"] ="<<rclk_a[i]<<endl;
      cout<<"d1c_a["<<i<<"] = "<<d1c_a[i]<<"   rclk["<<i<<"] ="<<rclk_a[i]<<endl;
    }
  
  //Define histos.
  TH1F *h1 = new TH1F("h1",Form("Beam Current U1 Run %i",run) , nevts+1, 0, nevts);
  TH1F *h2 = new TH1F("h2",Form("Beam Current D1 Run %i",run) , nevts+1, 0, nevts);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.08);
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetLabelSize(0.05);
  h1->GetYaxis()->SetTitleSize(0.06);
  h1->GetYaxis()->SetTitleOffset(0.75);
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetXaxis()->SetTitleSize(0.06);
  h1->GetXaxis()->SetTitleOffset(0.75);
  
  h2->GetYaxis()->CenterTitle(true);
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.06);
  h2->GetYaxis()->SetTitleOffset(0.75);
  h2->GetXaxis()->CenterTitle(true);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetXaxis()->SetTitleSize(0.06);
  h2->GetXaxis()->SetTitleOffset(0.75);
  
  // cout<<"!!!!!!!!!!!!!!!!"<<h1->GetBinNumber()<<endl;
  //Calculate and plot beam current over run period. BCM scalers update every 4 seconds.
  for (Int_t i=0; i<nevts; i++) 
    {
      //ev->GetEntry(i);
      //Iu1[i] = ( (u1c_a[i+1]-u1c_a[i])/(rclk_a[i+1]-rclk_a[i])-Koff_u1 ) / Cv2f_u1;
      //Id1[i] = ( (d1c_a[i+1]-d1c_a[i])/(rclk_a[i+1]-rclk_a[i])-Koff_d1 ) / Cv2f_d1;
      if(rclk_a[i+1]!=rclk_a[i])
	{
	  Iu1[i] = ( (u1c_a[i+1]-u1c_a[i])/((rclk_a[i+1]-rclk_a[i])/clk_r)-Koff_u1 ) / Cv2f_u1;
	  Id1[i] = ( (d1c_a[i+1]-d1c_a[i])/((rclk_a[i+1]-rclk_a[i])/clk_r)-Koff_d1 ) / Cv2f_d1;
	}
      else
	{
	  Iu1[i] = 0.;
	  Id1[i] = 0.;
	}

      //cout<<"Iu1["<<i<<"] = "<<Iu1[i]<<"   Id1["<<i<<"] = "<<Id1[i]<<endl;
      h1->SetBinContent(i+1,(Iu1[i]));
      h2->SetBinContent(i+1,(Id1[i]));
      //Could plot as 2D histo with i vs I[i] or as a TGraph which is probably best.
    }

  
  for (Int_t i=0; i<nevts; i++) 
    {
      //cout<<"Scaler Event # "<<i<<"   U1c = "<<u1c_a[i]<<"   rclk = "<<rclk_a[i]<<"   I = "<<I[i]<<" uA"<<endl;
    }

  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  c1->Divide(1,2);
  
  c1->cd(1);
  h1->GetXaxis()->SetTitle("Number of Events");
  h1->GetYaxis()->SetTitle("Beam Current (uA)");
  h1->Draw("L");

  //Now we need to remove beam trips from the charge calculation.
  Int_t thr_u1[nevts];    //Above or below set current minimum threshold.
  Double_t u1c_beam=0., rclk_beam_u1=0.;
  Int_t rising_edge_u1=0;
  
  for(Int_t i=0;i<nevts;i++)
    {
      if(Iu1[i]>=I_set)
	{
	  thr_u1[i] = 1;
	}
      else
	{
	  thr_u1[i] = 0;
	}
      cout<<"thr_u1["<<i<<"] = "<<thr_u1[i]<<endl;
    }

  for(Int_t i=1;i<nevts;i++)
    {
      //Could miss first step if thr=1 for i=0 and thr=0 for i=1. Add this bit to catch that possibility.
      if(thr_u1[0]==1 && thr_u1[1]==0)
	{
	  u1c_beam = u1c_beam + u1c_a[0];
	  rclk_beam_u1 = rclk_beam_u1 + rclk_a[0];
	}

      //Find falling edge of beam trip or final entry.
      if(thr_u1[i]==1 && (thr_u1[i+1]==0 || i+1==nevts) )
	{
	  //Look back to find last rising edge of beam trip.
	  for(Int_t j=1;j<nevts;j++)
	    {
	      if(thr_u1[i-j]==0 || (i-j)==0)
		{
		  rising_edge_u1 = i-j+1;
		  //cout<<"rising_edge_u1 = "<<rising_edge_u1<<endl;
		  break;
		}
	    }

	  //Sum the total scaler values for the BCMs and clock when there was good beam.
	  //u1c_beam = u1c_beam + u1c_a[i] - u1c_a[rising_edge_u1];
	  //rclk_beam_u1 = rclk_beam_u1 + rclk_a[i] - rclk_a[rising_edge_u1];

	  //Sum the total scaler values for the BCMs and clock when there was good beam.
	  //Must avoid checking above top of arrays so add if statement.
	  if(i!=nevts-1)
	    {
	      if((rising_edge_u1+edge_off)<nevts)
		{
		  if(rising_edge_u1!=1)   //Not looking at start of run with beam above threshold.
		    {
		      if((rising_edge_u1+edge_off)>(i+1))   //Make sure rising edge always before falling.
			{
			  rising_edge_u1 = i+1;
			  edge_off = 0;
			}
		      u1c_beam = u1c_beam + u1c_a[i+1] - u1c_a[rising_edge_u1+edge_off];	  
		      rclk_beam_u1 = rclk_beam_u1 + rclk_a[i+1] - rclk_a[rising_edge_u1+edge_off];
		      output1<<rising_edge_u1+edge_off<<"   "<<u1c_a[rising_edge_u1+edge_off]<<"   "<<i+1<<"   "<<u1c_a[i+1]<<endl;
		    }
		  else   //Looking at start of run with beam above threshold.
		    {
		      u1c_beam = u1c_beam + u1c_a[i+1] - u1c_a[0];	  
		      rclk_beam_u1 = rclk_beam_u1 + rclk_a[i+1] - rclk_a[0];
		      output1<<0<<"   "<<u1c_a[0]<<"   "<<i+1<<"   "<<u1c_a[i+1]<<endl;
		    }
		}
	      else   //Protect against going over max array entry.
		{
		  u1c_beam = u1c_beam + u1c_a[i+1] - u1c_a[nevts-1];	  
		  rclk_beam_u1 = rclk_beam_u1 + rclk_a[i+1] - rclk_a[nevts-1];
		  output1<<nevts-1<<"   "<<u1c_a[nevts-1]<<"   "<<i+1<<"   "<<u1c_a[i+1]<<endl;
		}
	    }
	  else
	    {
	      if(rising_edge_u1+edge_off>nevts) //Make sure not to overflow array.
		{
		  edge_off=0;
		}
	      u1c_beam = u1c_beam + u1c_a[i] - u1c_a[rising_edge_u1+edge_off];	  
	      rclk_beam_u1 = rclk_beam_u1 + rclk_a[i] - rclk_a[rising_edge_u1+edge_off]; 
	      output1<<rising_edge_u1+edge_off<<"   "<<u1c_a[rising_edge_u1+edge_off]<<"   "<<i<<"   "<<u1c_a[i]<<endl;
	    }

	  //cout<<"i = "<<i<<"   u1c_beam = "<<u1c_beam<<"   rclk_beam_u1 = "<<rclk_beam_u1<<endl;
	}

      //Draw beam trip cuts.
      if(thr_u1[i]!=thr_u1[i-1])
	{
	  //Separate cases for rising edge and falling edge. Only delay count on rising edge waiting for beam to stabilize.
	  c1->Update();
	  if(thr_u1[i]==0) //Falling edge
	    {
	      TLine *line = new TLine(i,gPad->GetUymin(),i,gPad->GetUymax());
	      line->SetLineColor(kRed);
	      line->Draw("same");
	    }
	  if(thr_u1[i]==1) //Rising edge
	    {
	      TLine *line = new TLine(i+edge_off,gPad->GetUymin(),i+edge_off,gPad->GetUymax());
	      line->SetLineColor(kRed);
	      line->Draw("same");
	    }
	}
      //cout<<"thr_u1["<<i<<"] = "<<thr_u1[i]<<endl;
    }

  cout<<"****************************************************************************************"<<endl;  

  //TLine *line = new TLine(50,0,50,h2->GetMaximum());
  //c1->Update();
  //TLine *line = new TLine(50,gPad->GetUymin(),50,gPad->GetUymax());
  //line->SetLineColor(kRed);
  //line->Draw();

  c1->cd(2);
  h2->GetXaxis()->SetTitle("Number of Events");
  h2->GetYaxis()->SetTitle("Beam Current (uA)");
  h2->Draw("L");

  //Now we need to remove beam trips from the charge calculation.
  Int_t thr_d1[nevts];    //Above or below set current minimum threshold.
  Double_t d1c_beam=0., rclk_beam_d1=0.;
  Int_t rising_edge_d1=0;
  
  for(Int_t i=0;i<nevts;i++)
    {
      if(Id1[i]>=I_set)
	{
	  thr_d1[i] = 1;
	}
      else
	{
	  thr_d1[i] = 0;
	}
      cout<<"thr_d1["<<i<<"] = "<<thr_d1[i]<<endl;
    }

  for(Int_t i=1;i<nevts;i++)
    {
      //Could miss first step if thr=1 for i=0 and thr=0 for i=1. Add this bit to catch that possibility.
      if(thr_d1[0]==1 && thr_d1[1]==0)
	{
	  d1c_beam = d1c_beam + d1c_a[0];
	  rclk_beam_d1 = rclk_beam_d1 + rclk_a[0];
	}

      //Find falling edge of beam trip or final entry.
      if(thr_d1[i]==1 && (thr_d1[i+1]==0 || i+1==nevts) )
	{
	  //Look back to find last rising edge of beam trip.
	  for(Int_t j=1;j<nevts;j++)
	    {
	      if(thr_d1[i-j]==0 || (i-j)==0)
		{
		  rising_edge_d1 = i-j+1;
		  cout<<"rising_edge_d1 = "<<rising_edge_d1<<endl;
		  break;
		}
	    }

	  //Sum the total scaler values for the BCMs and clock when there was good beam.
	  //Must avoid checking above top of arrays so add if statement.
	  if(i!=nevts-1)
	    {
	      if((rising_edge_d1+edge_off)<nevts)
		{
		  if(rising_edge_d1!=1)   //Not looking at start of run with beam above threshold.
		    {
		      if((rising_edge_d1+edge_off)>(i+1))   //Make sure rising edge always before falling.
			{
			  rising_edge_d1 = i+1;
			  edge_off = 0;
			}
		      d1c_beam = d1c_beam + d1c_a[i+1] - d1c_a[rising_edge_d1+edge_off];	  
		      rclk_beam_d1 = rclk_beam_d1 + rclk_a[i+1] - rclk_a[rising_edge_d1+edge_off];
		      output2<<rising_edge_d1+edge_off<<"   "<<d1c_a[rising_edge_d1+edge_off]<<"   "<<i+1<<"   "<<d1c_a[i+1]<<endl;
		    }
		  else   //Looking at start of run with beam above threshold.
		    {
		      d1c_beam = d1c_beam + d1c_a[i+1] - d1c_a[0];	  
		      rclk_beam_d1 = rclk_beam_d1 + rclk_a[i+1] - rclk_a[0];
		      output2<<0<<"   "<<d1c_a[0]<<"    "<<i+1<<"   "<<d1c_a[i+1]<<endl;
		    }
		}
	       else   //Protect against going over max array entry.
		{
		  d1c_beam = d1c_beam + d1c_a[i+1] - d1c_a[nevts-1];	  
		  rclk_beam_d1 = rclk_beam_d1 + rclk_a[i+1] - rclk_a[nevts-1];
		  output2<<nevts-1<<"   "<<d1c_a[nevts-1]<<"   "<<i+1<<"   "<<d1c_a[i+1]<<endl;
		}
	    }
	  else
	    {
	      d1c_beam = d1c_beam + d1c_a[i] - d1c_a[rising_edge_d1+edge_off];	  
	      rclk_beam_d1 = rclk_beam_d1 + rclk_a[i] - rclk_a[rising_edge_d1+edge_off];	
	      output2<<rising_edge_d1+edge_off<<"   "<<d1c_a[rising_edge_d1+edge_off]<<"   "<<i<<"   "<<d1c_a[i]<<endl;
	    }
	  //cout<<"i = "<<i<<"   d1c_beam = "<<d1c_beam<<"   rclk_beam_d1 = "<<rclk_beam_d1<<endl;
	}

      //Draw beam trip cuts.
      if(thr_d1[i]!=thr_d1[i-1])
	{

	  //Separate cases for rising edge and falling edge. Only delay count on rising edge waiting for beam to stabilize.
	  c1->Update();
	  if(thr_u1[i]==0) //Falling edge
	    {
	      TLine *line = new TLine(i,gPad->GetUymin(),i,gPad->GetUymax());
	      line->SetLineColor(kRed);
	      line->Draw("same");
	    }
	  if(thr_u1[i]==1) //Rising edge
	    {
	      TLine *line = new TLine(i+edge_off,gPad->GetUymin(),i+edge_off,gPad->GetUymax());
	      line->SetLineColor(kRed);
	      line->Draw("same");
	    }
	}
      //cout<<"thr_u1["<<i<<"] = "<<thr_u1[i]<<endl;
    }
  //cout<<thr_d1[232]<<"   nevts = "<<nevts<<endl;
  //Calculate charge taking into account the beam trips. 
  Double_t Qu1_tot, Qd1_tot;
  Double_t time_u1_tot, time_d1_tot;
  Double_t e_u1_tot, e_d1_tot;
  Double_t Q_avg, e_avg;

  time_u1_tot = rclk_beam_u1/clk_r;
  time_d1_tot = rclk_beam_d1/clk_r;
  Qu1_tot = ( (u1c_beam/time_u1_tot-Koff_u1)/(Cv2f_u1) )*time_u1_tot;
  Qd1_tot = ( (d1c_beam/time_d1_tot-Koff_d1)/(Cv2f_d1) )*time_d1_tot;
  
  e_u1_tot = Qu1_tot*pow(10,-6)/e;
  e_d1_tot = Qd1_tot*pow(10,-6)/e;

  Q_avg = (Qu1_tot + Qd1_tot)/2.;
  e_avg = (e_u1_tot + e_d1_tot)/2.;

  cout<<"Total Qu1 = "<<Qu1_tot<<" uC   Total Qd1 = "<<Qd1_tot<<" uC"<<endl;
  cout<<"Total e- U1 = "<<e_u1_tot<<"   Total e- D1 = "<<e_d1_tot<<endl;
  cout<<"Total charge average = "<<Q_avg<<" uC"<<"   Total e- average = "<<e_avg<<endl;
  output1.close();
  output2.close();
}
