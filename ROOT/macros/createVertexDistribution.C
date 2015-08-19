#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include <math.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <iostream>
#include <string>

#include <iostream>
#include <cmath>
#include "TLegend.h"
#include "TPad.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"

using namespace std;


void createVertexDistribution(){
  TString selection1="pt1>18&&pt2>20&&vertices>0&&diLeptons==0&&lPFIsoDB<0.1&&tauIso<2&&lTrigger==1&&tightMuons<=1";

  TFile *f    = new TFile("/nfs_scratch/laura/ztt_unweighted4/muDATA.root");
  TTree* tree = (TTree*)f->Get("muTauEventTree/eventTree");


  TH1F* Shape = new TH1F("Nominal","Nominal",25,0,50);
  tree->Draw("vertices>>+Nominal","("+selection1+")");
  Shape->Scale(1/Shape->Integral());

  TFile *f1    = new TFile("/nfs_scratch/laura/ztt_unweighted4/ZJets.root");
  TTree* tree1 = (TTree*)f1->Get("muTauEventTree/eventTree");
  TH1F* Shape1 = new TH1F("Nominal1","Nominal1",25,0,50);
  tree1->Draw("vertices>>+Nominal1","("+selection1+")");
  Shape1->Scale(1/Shape1->Integral());

  TFile f2("/nfs_scratch/laura/ztt_unweighted4/vertices.root","RECREATE");
  TH1F* Shape2 = new TH1F("Nominal2","Nominal2",25,0,50);
  Shape2->Divide(Shape,Shape1);
  std::cout<<"Bin content of bin with 12 vertices "<<Shape2->GetBinContent(Shape2->FindBin(20))<<std::endl;
  std::cout<<"Shape integral "<< Shape2->Integral()<<std::endl;
  //f2->cd();
  Shape2->Write("vertices",TObject::kOverwrite);//"vertices",TObject::kOverwrite);

}