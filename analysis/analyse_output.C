void analyse_output(const char * filename1, const char * filename2, const char * tag1, const char * tag2)
{
  TFile f1(filename1,"READ");
  TFile f2(filename2,"READ");

  vector<string> hists;
  hists.push_back("residual_vertices_x");
  hists.push_back("residual_vertices_y");
  hists.push_back("residual_vertices_z");
  hists.push_back("residual_vertices_t");
  hists.push_back("residual_vertices_phi");
  hists.push_back("residual_vertices_s");

  TH1F * h1 = 0, * h2 = 0;
  TCanvas c;
  c.SaveAs(TString::Format("comparision_%s_%s.pdf[", tag1, tag2));
  for(size_t ih = 0; ih < hists.size(); ih++) {
    f1.GetObject(hists[ih].c_str(), h1);
    f2.GetObject(hists[ih].c_str(), h2);

    cout << hists[ih] << endl
	 << tag1 << "\t" << h1->GetMean() << "\t+-\t" << h1->GetRMS() << endl
	 << tag2 << "\t" << h2->GetMean() << "\t+-\t" << h2->GetRMS() << endl
	 << endl;

    h1->SetLineColor(kRed);
    h2->SetLineStyle(kDashed);
    h1->Draw();
    h2->Draw("SAME");
    c.SaveAs(TString::Format("comparision_%s_%s.pdf", tag1, tag2));
    if(hists[ih] == "residual_vertices_t") {
      h1->GetXaxis()->SetRangeUser(-200,200);
      h1->Draw();
      h2->Draw("SAME");
      c.SaveAs(TString::Format("comparision_%s_%s.pdf", tag1, tag2));
    }
  }
  c.SaveAs(TString::Format("comparision_%s_%s.pdf]", tag1, tag2));
}
