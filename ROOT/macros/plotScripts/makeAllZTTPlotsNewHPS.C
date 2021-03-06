{

  gROOT->ProcessLine(".x UWAnalysis/ROOT/macros/plotters/loadMuTauPlotter.C");

  std::string lumi="35";

  //SELECTION REQUIREMENTS

  std::string selection = "((HLT_Mu9_wasRun==1&&HLT_Mu9_prescale==1&&HLT_Mu9_fired==1)||(HLT_Mu15_v1_wasRun==1&&HLT_Mu15_v1_prescale==1&&HLT_Mu15_v1_fired==1))&&PVs>0&&mumuSize==0&&muTauRelPFIso<0.05&muTauMuonVeto&&muTauEleVeto&&muTauCharge==0&&muTauMt1<40&&muTauLooseIso&&muTauPt2>20&&muTauPt1>15";

 TCanvas * mass =   plotter->makeStackedPlot("muTauMass",selection,lumi,20,0,200,"visible Mass","GeV/c^{2}",90,80,60,0.6,0.6,false,0.1,180,"mass",true);
 mass->SaveAs("mass.png");
 mass->SaveAs("mass.pdf");
 mass->SaveAs("mass.root");

}
