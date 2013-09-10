{
  gStyle->SetOptStat(0);
  
  file = TFile("weight_hist.root","READ");


  mr = hist_high_data;
  mr->SetTitle("Data In High R^{2}");
  mr_toy = hist_high_pred;
  mr_toy->SetTitle("Estimation In High R^{2}");

  
  //build the canvas
  TCanvas *c2 = new TCanvas("c2","example",1024,768);
  TPad *pad3 = new TPad("pad3","pad3",0,0.4,1,1);

  pad3->SetTopMargin(.2);
  pad3->SetBottomMargin(0);
  pad3->Draw();
  pad3->cd();

  pad3->SetLogy();
  mr->GetYaxis()->SetTitleSize(.12);
  mr->GetYaxis()->SetTitleOffset(.45);
  mr->GetYaxis()->SetLabelSize(.08);
  mr->GetYaxis()->SetNdivisions(8);
  mr->SetLineColor(kBlack);
  mr->SetMarkerSize(1.6);
  mr->SetMarkerStyle(20);
  copy = mr->DrawCopy("p");

  if(hist_signal->GetEntries() > 0){
    hist_signal->SetFillColor(kSpring);
    hist_signal->SetLineColor(kSpring);
    hist_signal->SetFillStyle(3001);
    hist_signal->Draw("same");
    copy->Draw("same");
  }


  mr_toy->SetLineWidth(0);
  mr_toy->SetLineColor(kAzure+10);
  mr_toy->SetFillColor(kAzure+10);
  mr_toy->SetFillStyle(3001);
  mr_toy->Draw("e2same");
  mr_toy->SetMarkerStyle(4);




  TPaveText *myText=new TPaveText(0.347,0.8063,.898,.9121,"NDC");
  //  TPaveText *myText = new TPaveText(0.2,0.7,0.4,0.);
  myText->SetTextFont(42.);
  myText->SetTextSize(0.06);
  myText->SetFillColor(0); 
  myText->SetTextAlign(12);
  myText->AddText("CMS Preliminary #sqrt{s} = 8 TeV #intL dt = 19.789 fb^{-1}");   
  myText->Draw("same");

  c2->cd();


  TPad *pad4 = new TPad("pad4","pad4",0,0,1,0.4);

  pad4->SetTopMargin(0);
  pad4->SetBottomMargin(.4);
  pad4->Draw();
  pad4->cd();

  mr->SetStats(0);

  mr->GetYaxis()->SetTitle("Data / Est.");
  mr->GetYaxis()->SetTitleOffset(.31);
  mr->GetYaxis()->SetTitleSize(.175);
  mr->GetYaxis()->SetLabelSize(.12);

  mr->GetXaxis()->SetLabelSize(.12);
  mr->GetXaxis()->SetTitleSize(.175);
  mr->GetXaxis()->SetTitle("M_{R} [TeV]");

  mr->Divide(mr_toy);
  mr->SetMarkerStyle(20);
  mr->SetFillColor(kOrange+1);
  mr->SetFillStyle(3001);  
  mr_toy->SetFillStyle(3001);
  mr->SetMaximum(1.5);
  mr->SetMinimum(.2);

  mr->Draw("e2p");
  line2 = TLine(mr->GetXaxis()->GetXmin(),1,mr->GetXaxis()->GetXmax(),1);
  line2.Draw("same");


  c2->cd();

  gStyle->SetOptTitle(0);
  leg3 = pad3->BuildLegend(.415,.394,.86,.738);
  //  leg3 = pad3->BuildLegend(1.406,2.22,2.89,3.941,"NDC");
  /*  leg3->SetX1(1.406);
  leg3->SetX2(2.89);
  leg3->SetY1(2.22);
  leg3->SetY2(3.941);*/
  leg3->SetFillColor(0);
  leg3->SetLineColor(0);

  TCanvas *c3 = new TCanvas("c3","example",1024,768);

  
  hist_difference->SetMarkerStyle(20);
  /*hist_difference->SetFillColor(kOrange+1);
  hist_difference->SetFillStyle(3001);*/
  hist_difference->Draw("e");


  hist_difference->SetTitle("Data - Prediction");
  hist_difference->GetYaxis()->SetTitle("N Data - N Predicted");
  hist_difference->GetXaxis()->SetTitle("M_{R} [TeV]");

  hist_difference->GetYaxis()->SetTitleOffset(.9);
  hist_difference->GetXaxis()->SetTitleOffset(.9);
  hist_difference->GetYaxis()->SetTitleSize(.075);
  hist_difference->GetYaxis()->SetLabelSize(.06);
  hist_difference->GetXaxis()->SetLabelSize(.06);
  hist_difference->GetXaxis()->SetTitleSize(.075);

  TPaveText *cms_2=new TPaveText(0.3,0.88,.881,.95,"NDC");

  //  TPaveText *cms_2 = new TPaveText(0.2,0.7,0.4,0.);
  cms_2->SetTextFont(42.);
  cms_2->SetTextSize(0.04);
  cms_2->SetFillColor(0); 
  cms_2->SetTextAlign(12);
  cms_2->AddText("CMS Preliminary #sqrt{s} = 8 TeV #intL dt = 19.789 fb^{-1}");   
  cms_2->Draw("same");


  

  c3->SetTopMargin(.12);

  if(hist_signal->GetEntries() > 0) {
    hist_signal->Draw("same");
    hist_difference->Draw("esame");
  }


  gStyle->SetOptTitle(0);

  if(hist_signal.GetEntries() > 0){
    leg4 = c3->BuildLegend();
    leg4->SetFillColor(0);
    leg4->SetLineColor(0);
  }

  line3 = TLine(hist_difference->GetXaxis()->GetXmin(),0,hist_difference->GetXaxis()->GetXmax(),0);
  line3.Draw("same");

  canvas_bin0->Draw();
}

   
