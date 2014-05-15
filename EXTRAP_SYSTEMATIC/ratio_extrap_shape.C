{
  //study_string_high = "Razor #gamma#gamma + #geq 1 jet: Signal Injection Study";
  //study_string_low = "Razor #gamma#gamma + #geq 1 jet: Signal Injection Study";
  
  //study_string_low = "Razor #gamma#gamma + #geq 1 jet: Low R^{2} Fit ff Sample";
  //study_string_high = "Razor #gamma#gamma + #geq 1 jet: High R^{2} ff Sample";
  study_string_low = "Razor #gamma#gamma + #geq 1 jet: Low R^{2} Fit #gamma#gamma Sample";  
  study_string_high = "Razor #gamma#gamma + #geq 1 jet: High R^{2} #gamma#gamma Sample";
  gStyle->SetOptStat(0);
  
  file = TFile("weight_hist.root","READ");

  mr = hist_high_data;
  //  mr->SetTitle("Data In High R^{2}");
  mr_toy = hist_high_pred;
  mr_toy->SetTitle("Estimation In High R^{2}");

  Float_t xmin = mr_toy->GetXaxis()->GetXmin();
  Float_t xmax = mr_toy->GetXaxis()->GetXmax();
  
  //build the canvas
  TCanvas *high_canvas = new TCanvas("c2","high canvas",1024,768);
  TPad *high_plot = new TPad("pad3","pad3",0,0.4,1,1);

  high_plot->SetTopMargin(.2);
  high_plot->SetBottomMargin(0);
  high_plot->Draw();
  high_plot->cd();

  high_plot->SetLogy();

  mr->SetLineColor(kBlack);
  mr->SetMarkerSize(1.2);
  mr->SetMarkerStyle(20);
  mr->GetYaxis()->SetTitleSize(.12);
  mr->GetYaxis()->SetTitleOffset(.45);
  mr->GetYaxis()->SetLabelSize(.13);
  mr->GetYaxis()->SetNdivisions(8);
  //copy = mr->DrawCopy("p");

  //  hist_high_up->SetFillStyle(3001);
  hist_high_up->SetFillColor(kAzure-9);
  hist_high_up->SetLineColor(kAzure-9);
  hist_high_up->SetLineWidth(2);
  //highgraph->Draw("a");
  //  

  //hist_high_down->Draw();


  //hist_high_up->Draw("same");

   

  //highgraph->Delete();


  //  mr_toy->SetLineWidth(2.0);
  //mr_toy->SetLineColor(kBlack);
  //mr_toy->SetFillColor(kAzure-9);
  //mr_toy->SetFillStyle(3000);

  //mr_toy->Draw("psame");
  mr_toy->SetMarkerStyle(4);
  fit_graph_high->GetXaxis()->SetRangeUser(xmin,xmax);
  fit_graph_high->GetYaxis()->SetLabelSize(.08);
  fit_graph_high->GetYaxis()->SetTitleSize(.09);
  fit_graph_high->GetYaxis()->SetTitleOffset(.6);
  fit_graph_high->GetYaxis()->SetTitle("N Events");
  fit_graph_high->SetLineColor(kAzure-9);
  //  fit_graph_high->SetFillStyle(3001);
  fit_graph_high->Draw("ep2a");
  mr_toy->Draw("psame");
    mr->Draw("e1psame"); 
    //hist_high_down->SetFillStyle(3000);
  hist_high_down->SetFillColor(kWhite-9);
  hist_high_down->SetLineColor(kAzure-9);
  //  hist_high_down->SetLineColor(kAzure-9);
  hist_high_down->SetLineWidth(2.0);
  //hist_high_down->Draw("same");
  //  hist_high_down->Draw("asame");
  


  high_plot->RedrawAxis();

  TPaveText *myText=new TPaveText(0.347,0.8063,.898,.9121,"NDC");
  //  TPaveText *myText = new TPaveText(0.2,0.7,0.4,0.);
  myText->SetTextFont(42.);
  myText->SetTextSize(0.07);
  myText->SetFillColor(0); 
  myText->SetTextAlign(12);
  myText->AddText("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}");   
  myText->Draw("same");


  Float_t xmin = mr_toy->GetXaxis()->GetXmin();
  Float_t xmax = mr_toy->GetXaxis()->GetXmax();

  mr->GetXaxis()->SetRangeUser(xmin,xmax);
  cout << xmin << " " <<  xmax << endl;

  gStyle->SetOptTitle(0);
  TLegend leg3(.58,.6,.93,.99,study_string_high,"NDC");
  leg3.SetTextSize(.07);
  leg3.AddEntry(fit_graph_high,"68% Shape Systematic","f");
  leg3.AddEntry(mr_toy,"Estimation in High R^{2}","p");
  leg3.AddEntry(mr,mr.GetTitle(),"pl");
  leg3.SetFillColor(0);
  leg3.SetLineWidth(3);


  if(hist_signal->GetEntries() > 0){
    //hist_signal->SetFillColor(kRed-9);
    hist_signal->SetLineColor(kRed-9);
    hist_signal->SetLineWidth(4);
    //hist_signal->SetFillStyle(3000);
    hist_signal->Draw("same");
    leg3.AddEntry(hist_signal,"Signal","l");
    mr->Draw("psame");
  }
  leg3.Draw("same");  
  
  high_canvas->cd();

  TPad *high_graph = new TPad("pad4","pad4",0,0,1,0.4);

  high_graph->SetTopMargin(0);
  high_graph->SetBottomMargin(.4);
  high_graph->Draw();
  high_graph->cd();

  mr->SetMarkerStyle(20);
  mr->Draw("e2p");
  
  //mr->Draw("");

  highgraph->Draw("e2pa");
  highgraph_obs->Draw("psame");
  highgraph->GetYaxis()->SetTitle("z-score");
  //highgraph->GetYaxis()->SetTitle("2(N_{obs}-N_{exp})/win_{68}");
  highgraph->GetYaxis()->SetTitleOffset(.42);
  highgraph->GetYaxis()->SetTitleSize(.13);
  highgraph->GetYaxis()->SetLabelSize(.13);
  highgraph->GetYaxis()->SetRangeUser(-2.5,2.5);
  highgraph->GetYaxis()->SetNdivisions(5);

  highgraph->GetXaxis()->SetLabelSize(.14);
  highgraph->GetXaxis()->SetTitleSize(.17);
  highgraph->GetXaxis()->SetTitle("M_{R} [TeV]");
  highgraph->GetXaxis()->SetLimits(xmin,xmax);
  highgraph->GetXaxis()->SetRangeUser(xmin,xmax);

  TLine * line3 = new TLine(xmin, 0, xmax, 0);
  line3->SetLineWidth(2);
  line3->SetLineStyle(2);

  line3->Draw("same");  
  
  high_canvas->cd();


   

  //  leg3 = high_plot->BuildLegend(.415,.394,.86,.738);
  //  leg3 = pad3->BuildLegend(1.406,2.22,2.89,3.941,"NDC");
  /*  leg3->SetX1(1.406);
  leg3->SetX2(2.89);
  leg3->SetY1(2.22);
  leg3->SetY2(3.941);*/
  leg3->SetFillColor(0);
  leg3->SetLineColor(0);

  TCanvas *c3 = new TCanvas("c3","difference",1024,768);

  
  hist_difference->SetMarkerStyle(20);
  /*hist_difference->SetFillColor(kOrange+1);
  hist_difference->SetFillStyle(3000);*/
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
  cms_2->AddText("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}");   
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

  TLine * line3  = new TLine(hist_difference->GetXaxis()->GetXmin(),0,hist_difference->GetXaxis()->GetXmax(),0);
  line3.Draw("same");

  high_canvas_bin0->Draw();


  mr_low = hist_low_data;
  //  mr_low->SetTitle("Data In Low R^{2}");
  mr_low_toy = hist_low_pred;
  mr_low_toy->SetTitle("Fit In Low R^{2}");
  // mr_low.SetLineWidth(0);
 // mr_low.Sumw2();
  //  mr_low_toy.Sumw2();
    
  
  //build the canvas
  TCanvas *low_canvas = new TCanvas("c4","low_canvas",1024,768);
  TPad *low_plot = new TPad("pad5","pad5",0,0.4,1,1);

  low_plot->SetTopMargin(.2);
  low_plot->SetBottomMargin(0);
  low_plot->Draw();
  low_plot->cd();

  low_plot->SetLogy();
  mr_low->GetYaxis()->SetTitleSize(.12);
  mr_low->GetYaxis()->SetTitleOffset(.45);
  mr_low->GetYaxis()->SetLabelSize(.08);
  mr_low->GetYaxis()->SetNdivisions(8);
  mr_low->SetLineColor(kBlack);
  mr_low->SetMarkerSize(1.2);
  mr_low->SetMarkerStyle(20);
  //  copy_low = mr_low->DrawCopy("p");


  //histogram for the low canvas
  mr = hist_low_data;
  //  mr->SetTitle("Data In Low R^{2}");
  mr_toy = hist_low_pred;
  mr_toy->SetTitle("Fit In Low R^{2}");

  //mr.Sumw2();
  //mr_toy.Sumw2();
    
  
  //build the low canvas
  /*  low_canvas->SetTopMargin(.2);
  low_canvas->SetBottomMargin(0);
  low_canvas->Draw();
  low_canvas->cd();
  */
  low_canvas->SetLogy();
  mr->GetYaxis()->SetTitleSize(.12);
  mr->GetYaxis()->SetTitleOffset(.45);
  mr->GetYaxis()->SetLabelSize(.08);
  mr->GetYaxis()->SetNdivisions(8);
  mr->SetLineColor(kBlack);
  //mr->SetLineWidth(0);
  mr->SetMarkerSize(1.2);
  mr->SetMarkerStyle(20);
  //  copy = mr->DrawCopy("p");

  //  mr_toy->SetLineWidth(0);
  //mr_toy->SetLineColor(kAzure-9);
  //mr_toy->SetFillColor(kAzure-9);
  //mr_toy->SetFillStyle(3000);
  
  hist_low_up->SetFillStyle(3000);
  hist_low_up->SetFillColor(kAzure-9);
  hist_low_up->SetLineColor(kAzure-9);
  hist_low_up->SetLineWidth(2);

  //hist_low_down->Draw();

  //hist_low_up->Draw("same");

  //mr_toy->Draw("psame");
  mr_toy->SetMarkerStyle(4);

  hist_low_down->SetFillStyle(3000);
  hist_low_down->SetFillColor(kBlack-9);
  hist_low_down->SetLineColor(kAzure-9);
  hist_low_down->SetLineWidth(2.0);
  //hist_low_down->Draw("same");

  fit_graph_low->GetXaxis()->SetRangeUser(xmin,xmax);
  fit_graph_low->GetYaxis()->SetLabelSize(.08);
  fit_graph_low->GetYaxis()->SetTitleSize(.09);
  fit_graph_low->GetYaxis()->SetTitleOffset(.6);
  fit_graph_low->GetYaxis()->SetTitle("N Events");
  fit_graph_low->SetLineColor(kAzure-9);
  //fit_graph_low->SetFillStyle(3000);
  fit_graph_low->Draw("ep2a");

  mr_toy->Draw("psame");
  mr_toy->SetMarkerStyle(4);

  low_plot->RedrawAxis();

  mr->Draw("e1psame");
  TPaveText *myText2=new TPaveText(0.347,0.8063,.898,.9121,"NDC");
  //  TPaveText *myText = new TPaveText(0.2,0.7,0.4,0.);
  myText2->SetTextFont(42.);
  myText2->SetTextSize(0.07);
  myText2->SetFillColor(0); 
  myText2->SetTextAlign(12);
  myText2->AddText("CMS Preliminary #sqrt{s} = 8 TeV, L  = 19.7 fb^{-1}");   
  myText2->Draw("same");

  TLegend leglow(.58,.6,.88,.95,study_string_low,"NDC");
  leglow.SetTextSize(.07);
  leglow.AddEntry(fit_graph_low,"68% Shape Systematic ","f");
  leglow.AddEntry(mr_toy,"Fit Prediction","p");
  leglow.AddEntry(mr,mr.GetTitle(),"pl");
  leglow.SetFillColor(0);
  leglow.SetLineWidth(0);
  leglow.SetLineColor(0);
  leglow.Draw("same");
  low_canvas->cd();

  TPad *low68 = new TPad("pad6","low68",0,0,1,0.4);

  low68->SetTopMargin(0);
  low68->SetBottomMargin(.4);
  low68->Draw();
  low68->cd();

  lowgraph->Draw("e2pa");
  lowgraph_obs->Draw("psame");
  lowgraph->GetYaxis()->SetTitle("z-score");

  //lowgraph->GetYaxis()->SetTitle("2(N_{obs}-N_{exp})/win_{68}");
  //  lowgraph->GetYaxis()->SetTitle("#Delta");
  lowgraph->GetYaxis()->SetTitleOffset(.42);
  lowgraph->GetYaxis()->SetTitleSize(.13);
  lowgraph->GetYaxis()->SetLabelSize(.13);
  lowgraph->GetYaxis()->SetRangeUser(-2.5,2.5);
  lowgraph->GetYaxis()->SetNdivisions(5);

  lowgraph->GetXaxis()->SetLabelSize(.14);
  lowgraph->GetXaxis()->SetTitleSize(.17);
  lowgraph->GetXaxis()->SetTitle("M_{R} [TeV]");
  lowgraph->GetXaxis()->SetLimits(xmin,xmax);
  lowgraph->GetXaxis()->SetRangeUser(xmin,xmax);

  TLine * line3 = new TLine(xmin, 0, xmax, 0);
  line3->SetLineWidth(2);
  line3->SetLineStyle(2);

  line3->Draw("same");  


  //  leg_low = low_plot->BuildLegend(.415,.394,.86,.738);
  //leg_low->SetFillColor(0);
  //leg_low->SetLineColor(0);
    
  TCanvas * ratio_canvas = new TCanvas("ratio","ratio",1024,768);

  cout << "attempting copy" << endl;

  hist_ratio->Draw();
  hist_ratio->GetXaxis()->SetRangeUser(-10,10);

  ratio_canvas->BuildLegend()->SetFillColor(0);

}

   
