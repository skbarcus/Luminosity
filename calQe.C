//=============Dien ===01/14/2017===============//

// This code uses to calculate the Qe for each //
//    target for each kin setting              //

//=======main key of the code =================//
//now I dont want to add branch in the tree 
//chain all the rootfile for one run and then
//use the bcm and clkcount information 
//to calculate Qe directly inside the code
//loop over the scaler event and then I calculate the 
//current for that event use that current as condition 
//to do the cut on beam trip cut.
//this way can be quick calculation.
//but to apply the beam trip cuts in analysis
// still need to add them to the tree

#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>

void calQe(){

  //set constant for calculation. This is from BCM calibration information 

  double u1c = 1035;
  double u1off = 200;
  double d1c = 1263;
  double d1off = 37;
  double clk_rate = 1024;
  const double beamtrip =10;

  std::string target; 
  std::string kin; 
  int run=0, lastrun=0;

  double P0, angle;
  std::string armL,armR;

  cerr<< " which target: "; cin >> target;
  cerr<< " which kin: ";  cin >> kin;

  //by this way this code can be use for any target and any kin
  //as long as it access to the same format of runlist file
  //so runlist file name need to carry Target and kin infor


  ifstream infile; infile.open(Form("../runlist_full/%s_%s_runlist.txt",target.data(),kin.data()));
  ofstream outfile ; outfile.open(Form("../%s/%s_%s_Qe_full.dat",kin.data(), target.data(), kin.data()));

  // const TString rootfilePath = "/work/halla/e08014/disk1/Rootfiles_Jan092017/";
  // const TString rootfilePath = "/work/halla/e08014/disk1/Rootfiles/";
  const TString rootfilePath =  "/work/halla/e08014/disk1/dien/E08014/replay_me/Rootfiles_official/";


  while(true){

    infile >> run>> kin >> P0 >> angle >> target>> armL >> armR;
    cout<< " check the file: "<< infile.eof()<<endl;

    if(infile.eof())break;



    cerr << run << "   "<< kin<< "   "<< P0<< "   "<< angle<< "  "<< target<< "  "<< armL<< "   "<< armR<< endl;
  
    TChain *RIGHT = new TChain("RIGHT");
    
    //==to add the splits rootfiles to one chain=====//

    std::ostringstream str;
    str << rootfilePath<<"e08014_"<<run;
    TString basename = str.str().c_str();
    TString rootfile = basename + ".root";

   cout<< " check the rootfile : "<< rootfile<<endl;


    TFile *file = new TFile(rootfile,"r");
    if(file->IsZombie())){
    cout<<" this rootfile doest not exist: "<<endl;
  }



    Long_t i=0;
    while ( !gSystem->AccessPathName(rootfile.Data()) ) {
        RIGHT->Add(rootfile.Data());

        cout << "ROOT file " << rootfile << " added to TChain." << endl;
        i++;
        //rootfile = basename + Form("_%d",i) + ".root";
        rootfile = basename + "_" + i + ".root";

    }
    //=====End adding up the slpits rootfiles=======//

    double C_u1, C_d1,right_current,bcm_u1c,bcm_d1c;
    double clkcount;

    double time; 
    double u1_count;
    double d1_count;
    double clk_count;

    //set branch address for bcm branch 

   RIGHT->SetBranchAddress("bcm_u1c",&bcm_u1c);
   RIGHT->SetBranchAddress("bcm_d1c", &bcm_d1c);
   //  RIGHT->SetBranchAddress("C_u1", &C_u1);
   // RIGHT->SetBranchAddress("C_d1", &C_d1);
   //RIGHT->SetBranchAddress("right_current",&right_current);
   RIGHT->SetBranchAddress("clkcount",&clkcount);

  int Nentry=RIGHT->GetEntries();
  int N= Nentry;
  cout<<" n : "<< N<<endl;
  double sum_u=0;
  double sum_d=0;
  double count=0;
  double count_av=0;
  double sum_u_av=0;
  double sum_d_av=0;
  double Q_tot,Qe, I, I_avg;
  double sum_I=0;
  double prev_u1=0;
  double prev_clk=0;
  double prev_d1=0;
  double I_u,I_d;

for (int j=0; j<N;j++){

   RIGHT->GetEntry(j); 

   //this one to loop over and check the scaler event 
   //this to set condition for cut after calculation
   //we just calculate with event got updated in time 

   if (clkcount!=prev_clk){
    
     count++;
     u1_count =bcm_u1c-prev_u1;
     d1_count = bcm_d1c-prev_d1;
     time = (clkcount-prev_clk)/clk_rate;

     //calculate the current 

     I_u = (u1_count/time -u1off)/u1c;
     I_d = (d1_count/time -d1off)/d1c;
     right_current = (I_u+I_d)/2;

     cout<< " I_u, I_d: "<< I_u << "   "<< I_d<< endl;

     //calculate the charge

     C_u1 = I_u*time;
     C_d1 = I_d*time;

     //calculate the charge with beam trip cuts now

    if (right_current>beamtrip){

      sum_u=sum_u+C_u1;
      sum_d=sum_d+C_d1;

       }

     prev_u1=bcm_u1c; prev_d1= bcm_d1c;
     prev_clk = clkcount;
    
    }
    
 }
   Qe = (sum_u+sum_d)/2;

   outfile << setiosflags(ios::left) << setw(5) << run << "     ";

   outfile << setiosflags(ios::left) << setw(8) << Qe << endl;

  }



}