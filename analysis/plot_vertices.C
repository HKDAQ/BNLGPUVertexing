#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <TH2F.h>
#include <TFile.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TStyle.h>
//#include <T.h>

struct point
{
  int ev;
  double x;
  double y;
  double z;
  double t;
};  

#define MAXARRAY 1000000

void do_fit(double distance, TH1F * h);

int plot_vertices(const char * fitresultfile, const char * truevertexfile, const char * outtag){

  cout << "STARTING>.."<< endl;
  double detector_H = 6000.;
  double detector_R = 3710.46;
  double distance = 500.;
  int nR = 2*2*detector_R/distance;
  int nH = 2*detector_H/distance;

  TFile output(TString::Format("vertices_%s.root",outtag),"RECREATE");

  ifstream verticesfile;
  verticesfile.open(fitresultfile);

  int nevent;
  double vx, vy, vz, time;

  point vtx[MAXARRAY];
  int nvertex = 0;

  while(!verticesfile.eof()){
    verticesfile >> nevent;
    verticesfile >> vx;
    verticesfile >> vy;
    verticesfile >> vz;
    verticesfile >> time;

    point v;
    v.ev = nevent;
    v.x = vx;
    v.y = vy;
    v.z = vz;
    v.t = time;

    vtx[nvertex] = v;
    nvertex ++;

  } 
  verticesfile.close(); 


  clog << " n vertices: " << nvertex << " nR " << nR << " nH " << nH << endl;

  TH1F vertices_x("vertices_x", "vertices x; x", nR,1,-1);
  vertices_x.SetLineColor(kBlack);
  vertices_x.SetLineWidth(2);

  TH1F vertices_y("vertices_y", "vertices y; y", nR,1,-1);
  vertices_y.SetLineColor(kBlack);
  vertices_y.SetLineWidth(2);

  TH1F vertices_phi("vertices_phi", "vertices phi; phi", nR*3,1,-1);
  vertices_phi.SetLineColor(kBlack);
  vertices_phi.SetLineWidth(2);

  TH1F vertices_z("vertices_z", "vertices z; z", nH,1,-1);
  vertices_z.SetLineColor(kBlack);
  vertices_z.SetLineWidth(2);

  TH1F vertices_s("vertices_s", "vertices s; s", nR,1,-1);
  vertices_s.SetLineColor(kBlack);
  vertices_s.SetLineWidth(2);

  TH1F vertices_t("vertices_t", "vertices t; t", 100,1,-1);
  vertices_t.SetLineColor(kBlack);
  vertices_t.SetLineWidth(2);

  TH2F vertices_y_vs_x("vertices_y_vs_x", "vertices y vs x; x; y", 
		       nR,1,-1, nR,1,-1);

  TH2F vertices_z_vs_phi("vertices_z_vs_phi", "vertices z vs phi; #phi; z", 
		       nR*3,1,-1, nH,1,-1);

  for( int i=0; i<nvertex; i++){
    vertices_x.Fill(vtx[i].x);
    vertices_y.Fill(vtx[i].y);
    vertices_z.Fill(vtx[i].z);
    vertices_t.Fill(vtx[i].t);
    vertices_phi.Fill(180./acos(-1.)*atan2(vtx[i].y, vtx[i].x));
    vertices_s.Fill(sqrt(pow(vtx[i].x,2) + pow(vtx[i].y,2)));
    vertices_y_vs_x.Fill(vtx[i].x, vtx[i].y); 
    vertices_z_vs_phi.Fill(180./acos(-1.)*atan2(vtx[i].y, vtx[i].x), vtx[i].z);
  }



  ifstream trueverticesfile;
  trueverticesfile.open(truevertexfile);

  point true_vtx[MAXARRAY];
  int ntruevertex = 0;

  while(!trueverticesfile.eof()){
    trueverticesfile >> vx;
    trueverticesfile >> vy;
    trueverticesfile >> vz;

    point v;
    v.ev = ntruevertex;
    v.x = vx;
    v.y = vy;
    v.z = vz;
    v.t = 0.;

    true_vtx[ntruevertex] = v;
    ntruevertex++;

  } 
  trueverticesfile.close(); 


  clog << " n true vertices: " << ntruevertex << " nR " << nR << " nH " << nH << endl;

  TH1F true_vertices_x("true_vertices_x", "true_vertices x; x", nR,1,-1);
  true_vertices_x.SetLineColor(kBlack);
  true_vertices_x.SetLineWidth(2);

  TH1F true_vertices_y("true_vertices_y", "true_vertices y; y", nR,1,-1);
  true_vertices_y.SetLineColor(kBlack);
  true_vertices_y.SetLineWidth(2);

  TH1F true_vertices_phi("true_vertices_phi", "true_vertices phi; phi", nR*3,1,-1);
  true_vertices_phi.SetLineColor(kBlack);
  true_vertices_phi.SetLineWidth(2);

  TH1F true_vertices_z("true_vertices_z", "true_vertices z; z", nH,1,-1);
  true_vertices_z.SetLineColor(kBlack);
  true_vertices_z.SetLineWidth(2);

  TH1F true_vertices_s("true_vertices_s", "true_vertices s; s", nR,1,-1);
  true_vertices_s.SetLineColor(kBlack);
  true_vertices_s.SetLineWidth(2);

  TH2F true_vertices_y_vs_x("true_vertices_y_vs_x", "true_vertices y vs x; x; y", 
		       nR,1,-1, nR,1,-1);

  TH2F true_vertices_z_vs_phi("true_vertices_z_vs_phi", "true_vertices z vs phi; #phi; z", 
		       nR*3,1,-1, nH,1,-1);


  TH1F residual_vertices_x("residual_vertices_x", "residual_vertices x; x", 10*nR,1,-1);
  residual_vertices_x.SetLineColor(kBlack);
  residual_vertices_x.SetLineWidth(2);

  TH1F * residual_vertices_x_zone[9];
  TH1F * residual_vertices_y_zone[9];
  TH1F * residual_vertices_t_zone[9];

  residual_vertices_x_zone[0] = new TH1F("residual_vertices_x_xp_yp", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_x_zone[1] = new TH1F("residual_vertices_x_xp_y0", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_x_zone[2] = new TH1F("residual_vertices_x_xp_yn", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_x_zone[3] = new TH1F("residual_vertices_x_x0_yp", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_x_zone[4] = new TH1F("residual_vertices_x_x0_y0", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_x_zone[5] = new TH1F("residual_vertices_x_x0_yn", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_x_zone[6] = new TH1F("residual_vertices_x_xn_yp", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_x_zone[7] = new TH1F("residual_vertices_x_xn_y0", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_x_zone[8] = new TH1F("residual_vertices_x_xn_yn", "residual_vertices x; residual x", 10*nR,1,-1);
  residual_vertices_y_zone[0] = new TH1F("residual_vertices_y_xp_yp", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_y_zone[1] = new TH1F("residual_vertices_y_xp_y0", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_y_zone[2] = new TH1F("residual_vertices_y_xp_yn", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_y_zone[3] = new TH1F("residual_vertices_y_x0_yp", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_y_zone[4] = new TH1F("residual_vertices_y_x0_y0", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_y_zone[5] = new TH1F("residual_vertices_y_x0_yn", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_y_zone[6] = new TH1F("residual_vertices_y_xn_yp", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_y_zone[7] = new TH1F("residual_vertices_y_xn_y0", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_y_zone[8] = new TH1F("residual_vertices_y_xn_yn", "residual_vertices y; residual y", 10*nR,1,-1);
  residual_vertices_t_zone[0] = new TH1F("residual_vertices_t_xp_yp", "residual_vertices t; residual t", 10*nR,1,-1);
  residual_vertices_t_zone[1] = new TH1F("residual_vertices_t_xp_y0", "residual_vertices t; residual t", 10*nR,1,-1);
  residual_vertices_t_zone[2] = new TH1F("residual_vertices_t_xp_yn", "residual_vertices t; residual t", 10*nR,1,-1);
  residual_vertices_t_zone[3] = new TH1F("residual_vertices_t_x0_yp", "residual_vertices t; residual t", 10*nR,1,-1);
  residual_vertices_t_zone[4] = new TH1F("residual_vertices_t_x0_y0", "residual_vertices t; residual t", 10*nR,1,-1);
  residual_vertices_t_zone[5] = new TH1F("residual_vertices_t_x0_yn", "residual_vertices t; residual t", 10*nR,1,-1);
  residual_vertices_t_zone[6] = new TH1F("residual_vertices_t_xn_yp", "residual_vertices t; residual t", 10*nR,1,-1);
  residual_vertices_t_zone[7] = new TH1F("residual_vertices_t_xn_y0", "residual_vertices t; residual t", 10*nR,1,-1);
  residual_vertices_t_zone[8] = new TH1F("residual_vertices_t_xn_yn", "residual_vertices t; residual t", 10*nR,1,-1);
  for(int i=0; i<9; i++){
    residual_vertices_x_zone[i]->SetLineColor(kBlack);
    residual_vertices_x_zone[i]->SetLineWidth(2);
    residual_vertices_y_zone[i]->SetLineColor(kBlack);
    residual_vertices_y_zone[i]->SetLineWidth(2);
    residual_vertices_t_zone[i]->SetLineColor(kBlack);
    residual_vertices_t_zone[i]->SetLineWidth(2);
  }


  TH1F residual_vertices_y("residual_vertices_y", "residual_vertices y; y", 10*nR,1,-1);
  residual_vertices_y.SetLineColor(kBlack);
  residual_vertices_y.SetLineWidth(2);

  TH1F residual_vertices_y_pos_x("residual_vertices_y_pos_x", "residual_vertices y; y", 10*nR,1,-1);
  residual_vertices_y_pos_x.SetLineColor(kBlack);
  residual_vertices_y_pos_x.SetLineWidth(2);

  TH1F residual_vertices_y_neg_x("residual_vertices_y_neg_x", "residual_vertices y; y", 10*nR,1,-1);
  residual_vertices_y_neg_x.SetLineColor(kBlack);
  residual_vertices_y_neg_x.SetLineWidth(2);

  TH1F residual_vertices_phi("residual_vertices_phi", "residual_vertices phi; phi", 10*nR*3,1,-1);
  residual_vertices_phi.SetLineColor(kBlack);
  residual_vertices_phi.SetLineWidth(2);

  TH1F residual_vertices_z("residual_vertices_z", "residual_vertices z; z", 10*nH,1,-1);
  residual_vertices_z.SetLineColor(kBlack);
  residual_vertices_z.SetLineWidth(2);

  TH1F residual_vertices_t("residual_vertices_t", "residual_vertices t; t", 1000,1,-1);
  residual_vertices_t.SetLineColor(kBlack);
  residual_vertices_t.SetLineWidth(2);

  TH1F residual_vertices_s("residual_vertices_s", "residual_vertices s; s", 10*nR,1,-1);
  residual_vertices_s.SetLineColor(kBlack);
  residual_vertices_s.SetLineWidth(2);

  TH1F residual_vertices_s_pos_x("residual_vertices_s_pos_x", "residual_vertices s; s", 10*nR,1,-1);
  residual_vertices_s_pos_x.SetLineColor(kBlack);
  residual_vertices_s_pos_x.SetLineWidth(2);

  TH1F residual_vertices_s_neg_x("residual_vertices_s_neg_x", "residual_vertices s; s", 10*nR,1,-1);
  residual_vertices_s_neg_x.SetLineColor(kBlack);
  residual_vertices_s_neg_x.SetLineWidth(2);

  TH2F residual_vertices_s__vs__s("residual_vertices_s__vs__s", "residual_vertices s__vs__s; s; residual s", nR,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_s_t("residual_vertices_s_t", "residual_vertices s_t; t; s", 1000,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_s_t_pos_x("residual_vertices_s_t_pos_x", "residual_vertices s_t; t; s", 1000,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_s_t_neg_x("residual_vertices_s_t_neg_x", "residual_vertices s_t; t; s", 1000,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_s__vs__z("residual_vertices_s__vs__z", "residual_vertices s__vs__z; residual z; residual s", 10*nH,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_s__vs__x("residual_vertices_s__vs__x", "residual_vertices s__vs__x; x; residual s", nR,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_x__vs__x("residual_vertices_x__vs__x", "residual_vertices x__vs__x; x; residual x", nR,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_y__vs__x("residual_vertices_y__vs__x", "residual_vertices y__vs__x; x; residual y", nR,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_y__vs__y_pos_x("residual_vertices_y__vs__y_pos_x", "residual_vertices y__vs__y; y; residual y", nR,1,-1, 10*nR,1,-1);

  TH2F residual_vertices_y__vs__y_neg_x("residual_vertices_y__vs__y_neg_x", "residual_vertices y__vs__x; x; residual y", nR,1,-1, 10*nR,1,-1);


  double limit = 3700./4.;

  int ntruthmatched = 0;
  
  for( int i=0; i<ntruevertex; i++){
    true_vertices_x.Fill(true_vtx[i].x);
    true_vertices_y.Fill(true_vtx[i].y);
    true_vertices_z.Fill(true_vtx[i].z);
    true_vertices_phi.Fill(180./acos(-1.)*atan2(true_vtx[i].y, true_vtx[i].x));
    true_vertices_s.Fill(sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
    true_vertices_y_vs_x.Fill(true_vtx[i].x, true_vtx[i].y); 
    true_vertices_z_vs_phi.Fill(180./acos(-1.)*atan2(true_vtx[i].y, true_vtx[i].x), true_vtx[i].z);

    int zone = -1;

    if( true_vtx[i].x > limit && true_vtx[i].y > limit ){
      zone = 0;
    }else if( true_vtx[i].x > limit && fabs(true_vtx[i].y) < limit ){
      zone = 1;
    }else if( true_vtx[i].x > limit && true_vtx[i].y < -limit ){
      zone = 2;
    }else if( fabs(true_vtx[i].x) < limit && true_vtx[i].y > limit ){
      zone = 3;
    }else if( fabs(true_vtx[i].x) < limit && fabs(true_vtx[i].y) < limit ){
      zone = 4;
    }else if( fabs(true_vtx[i].x) < limit && true_vtx[i].y < -limit ){
      zone = 5;
    }else if( true_vtx[i].x < -limit && true_vtx[i].y > limit ){
      zone = 6;
    }else if( true_vtx[i].x < -limit && fabs(true_vtx[i].y) < limit ){
      zone = 7;
    }else if( true_vtx[i].x < -limit && true_vtx[i].y < -limit ){
      zone = 8;
    }
    
  
    for( int j=0; j<nvertex; j++){
      
      if( true_vtx[i].ev == vtx[j].ev ){
	//do some crude truth matching
	double distancediff = sqrt(pow(vtx[j].x - true_vtx[i].x, 2) + pow(vtx[j].y - true_vtx[i].y, 2) + pow(vtx[j].z - true_vtx[i].z, 2));
	double timediff = abs(vtx[j].t - true_vtx[i].t);
	if(timediff < 500 && distancediff < 1500)
	  ntruthmatched++;

	residual_vertices_x.Fill(vtx[j].x - true_vtx[i].x);

	residual_vertices_x_zone[zone]->Fill(vtx[j].x - true_vtx[i].x);
	residual_vertices_y_zone[zone]->Fill(vtx[j].y - true_vtx[i].y);
	residual_vertices_t_zone[zone]->Fill(vtx[j].t - true_vtx[i].t);

	residual_vertices_y.Fill(vtx[j].y - true_vtx[i].y);
	if( true_vtx[i].x > 0. )
	  residual_vertices_y_pos_x.Fill(vtx[j].y - true_vtx[i].y);
	else
	  residual_vertices_y_neg_x.Fill(vtx[j].y - true_vtx[i].y);
	residual_vertices_z.Fill(vtx[j].z - true_vtx[i].z);
	residual_vertices_t.Fill(vtx[j].t - true_vtx[i].t);
	residual_vertices_s.Fill(sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	if( true_vtx[i].x > 0. )
	  residual_vertices_s_pos_x.Fill(sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	else
	  residual_vertices_s_neg_x.Fill(sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	residual_vertices_s__vs__s.Fill( sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)),
					 sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	residual_vertices_s__vs__x.Fill( true_vtx[i].x,
					 sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	residual_vertices_x__vs__x.Fill( true_vtx[i].x,
					 vtx[j].x - true_vtx[i].x);
	residual_vertices_y__vs__x.Fill( true_vtx[i].x,
					 vtx[j].y - true_vtx[i].y);
	if( true_vtx[i].x > 0. )
	  residual_vertices_y__vs__y_pos_x.Fill( true_vtx[i].y,
					   vtx[j].y - true_vtx[i].y);
	else
	  residual_vertices_y__vs__y_neg_x.Fill( true_vtx[i].y,
					   vtx[j].y - true_vtx[i].y);
	residual_vertices_s_t.Fill(vtx[j].t - true_vtx[i].t,
				   sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	if( true_vtx[i].x > 0. )
	  residual_vertices_s_t_pos_x.Fill(vtx[j].t - true_vtx[i].t,
				     sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	else
	  residual_vertices_s_t_neg_x.Fill(vtx[j].t - true_vtx[i].t,
				     sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	residual_vertices_s__vs__z.Fill(vtx[j].z - true_vtx[i].z
					,
					sqrt(pow(vtx[j].x,2) + pow(vtx[j].y,2)) - sqrt(pow(true_vtx[i].x,2) + pow(true_vtx[i].y,2)));
	residual_vertices_phi.Fill(180./acos(-1.)*atan2(vtx[j].y, vtx[j].x) - 180./acos(-1.)*atan2(true_vtx[i].y, true_vtx[i].x));
      }
      if( true_vtx[i].ev < vtx[j].ev ) break;
    }

  }

  cout << "TRUTH MATCHING\t" << ntruthmatched << endl;

  TProfile* p_residual_vertices_s_t = residual_vertices_s_t.ProfileX(); 
  p_residual_vertices_s_t->GetXaxis()->SetTitle(" residual t [ns]");
  p_residual_vertices_s_t->GetYaxis()->SetTitle(" residual s [cm]");

  TProfile* p_residual_vertices_s_t_pos_x = residual_vertices_s_t_pos_x.ProfileX(); 
  p_residual_vertices_s_t_pos_x->GetXaxis()->SetTitle(" residual t [ns]");
  p_residual_vertices_s_t_pos_x->GetYaxis()->SetTitle(" residual s [cm]");

  TProfile* p_residual_vertices_s_t_neg_x = residual_vertices_s_t_neg_x.ProfileX(); 
  p_residual_vertices_s_t_neg_x->GetXaxis()->SetTitle(" residual t [ns]");
  p_residual_vertices_s_t_neg_x->GetYaxis()->SetTitle(" residual s [cm]");

  TProfile* p_residual_vertices_s__vs__s = residual_vertices_s__vs__s.ProfileX(); 
  p_residual_vertices_s__vs__s->GetXaxis()->SetTitle(" s [cm]");
  p_residual_vertices_s__vs__s->GetYaxis()->SetTitle(" residual s [cm]");

  TProfile* p_residual_vertices_s__vs__x = residual_vertices_s__vs__x.ProfileX(); 
  p_residual_vertices_s__vs__x->GetXaxis()->SetTitle(" true x [cm]");
  p_residual_vertices_s__vs__x->GetYaxis()->SetTitle(" residual s [cm]");

  TProfile* p_residual_vertices_x__vs__x = residual_vertices_x__vs__x.ProfileX(); 
  p_residual_vertices_x__vs__x->GetXaxis()->SetTitle(" true x [cm]");
  p_residual_vertices_x__vs__x->GetYaxis()->SetTitle(" residual x [cm]");

  TProfile* p_residual_vertices_y__vs__x = residual_vertices_y__vs__x.ProfileX(); 
  p_residual_vertices_y__vs__x->GetXaxis()->SetTitle(" true x [cm]");
  p_residual_vertices_y__vs__x->GetYaxis()->SetTitle(" residual y [cm]");

  TProfile* p_residual_vertices_y__vs__y_pos_x = residual_vertices_y__vs__y_pos_x.ProfileX(); 
  p_residual_vertices_y__vs__y_pos_x->GetXaxis()->SetTitle(" true y [cm]");
  p_residual_vertices_y__vs__y_pos_x->GetYaxis()->SetTitle(" residual y [cm]");

  TProfile* p_residual_vertices_y__vs__y_neg_x = residual_vertices_y__vs__y_neg_x.ProfileX(); 
  p_residual_vertices_y__vs__y_neg_x->GetXaxis()->SetTitle(" true y [cm]");
  p_residual_vertices_y__vs__y_neg_x->GetYaxis()->SetTitle(" residual y [cm]");

  TH1F * ratio_vertices_x = (TH1F*)vertices_x.Clone("ratio_vertices_x");
  ratio_vertices_x->Divide(&true_vertices_x);

  TH1F * ratio_vertices_y = (TH1F*)vertices_y.Clone("ratio_vertices_y");
  ratio_vertices_y->Divide(&true_vertices_y);

  TH1F * ratio_vertices_z = (TH1F*)vertices_z.Clone("ratio_vertices_z");
  ratio_vertices_z->Divide(&true_vertices_z);

  TH1F * ratio_vertices_phi = (TH1F*)vertices_phi.Clone("ratio_vertices_phi");
  ratio_vertices_phi->Divide(&true_vertices_phi);

  TH1F * ratio_vertices_s = (TH1F*)vertices_s.Clone("ratio_vertices_s");
  ratio_vertices_s->Divide(&true_vertices_s);

  TH2F * ratio_vertices_y_vs_x = (TH2F*)vertices_y_vs_x.Clone("ratio_vertices_y_vs_x");
  ratio_vertices_y_vs_x->Divide(&true_vertices_y_vs_x);

  TH2F * ratio_vertices_z_vs_phi = (TH2F*)vertices_z_vs_phi.Clone("ratio_vertices_z_vs_phi");
  ratio_vertices_z_vs_phi->Divide(&true_vertices_z_vs_phi);

  vertices_x.Sumw2();
  vertices_x.Scale(1./vertices_x.GetMaximum());
  vertices_y.Sumw2();
  vertices_y.Scale(1./vertices_y.GetMaximum());
  vertices_z.Sumw2();
  vertices_z.Scale(1./vertices_z.GetMaximum());
  vertices_phi.Sumw2();
  vertices_phi.Scale(1./vertices_phi.GetMaximum());
  vertices_s.Sumw2();
  vertices_s.Scale(1./vertices_s.GetMaximum());
  vertices_y_vs_x.Sumw2();
  vertices_y_vs_x.Scale(1./vertices_y_vs_x.GetMaximum());
  vertices_z_vs_phi.Sumw2();
  vertices_z_vs_phi.Scale(1./vertices_z_vs_phi.GetMaximum());

  p_residual_vertices_s_t->Write();
  p_residual_vertices_s_t_pos_x->Write();
  p_residual_vertices_s_t_neg_x->Write();
  p_residual_vertices_s__vs__s->Write();
  p_residual_vertices_s__vs__x->Write();
  p_residual_vertices_x__vs__x->Write();
  p_residual_vertices_y__vs__x->Write();
  p_residual_vertices_y__vs__y_pos_x->Write();
  p_residual_vertices_y__vs__y_neg_x->Write();

  output.Write();
  output.Close();

  TFile input(TString::Format("vertices_%s.root",outtag),"READ");

  TH1F * _vertices_x = (TH1F*)input.Get("vertices_x");
  TH1F * _vertices_y = (TH1F*)input.Get("vertices_y");
  TH1F * _vertices_z = (TH1F*)input.Get("vertices_z");
  TH1F * _vertices_phi = (TH1F*)input.Get("vertices_phi");
  TH1F * _vertices_s = (TH1F*)input.Get("vertices_s");
  TH2F * _vertices_y_vs_x = (TH2F*)input.Get("vertices_y_vs_x");
  TH2F * _vertices_z_vs_phi = (TH2F*)input.Get("vertices_z_vs_phi");

  TH1F * _true_vertices_x = (TH1F*)input.Get("true_vertices_x");
  TH1F * _true_vertices_y = (TH1F*)input.Get("true_vertices_y");
  TH1F * _true_vertices_z = (TH1F*)input.Get("true_vertices_z");
  TH1F * _true_vertices_phi = (TH1F*)input.Get("true_vertices_phi");
  TH1F * _true_vertices_s = (TH1F*)input.Get("true_vertices_s");
  TH2F * _true_vertices_y_vs_x = (TH2F*)input.Get("true_vertices_y_vs_x");
  TH2F * _true_vertices_z_vs_phi = (TH2F*)input.Get("true_vertices_z_vs_phi");

  TH1F * _ratio_vertices_x = (TH1F*)input.Get("ratio_vertices_x");
  TH1F * _ratio_vertices_y = (TH1F*)input.Get("ratio_vertices_y");
  TH1F * _ratio_vertices_z = (TH1F*)input.Get("ratio_vertices_z");
  TH1F * _ratio_vertices_phi = (TH1F*)input.Get("ratio_vertices_phi");
  TH1F * _ratio_vertices_s = (TH1F*)input.Get("ratio_vertices_s");
  TH2F * _ratio_vertices_y_vs_x = (TH2F*)input.Get("ratio_vertices_y_vs_x");
  TH2F * _ratio_vertices_z_vs_phi = (TH2F*)input.Get("ratio_vertices_z_vs_phi");

  TH2F * _residual_vertices_s_t = (TH2F*)input.Get("residual_vertices_s_t");
  TH2F * _residual_vertices_s__vs__s = (TH2F*)input.Get("residual_vertices_s__vs__s");
  TH2F * _residual_vertices_s__vs__z = (TH2F*)input.Get("residual_vertices_s__vs__z");

  TH1F * _residual_vertices_x = (TH1F*)input.Get("residual_vertices_x");
  TH1F * _residual_vertices_y = (TH1F*)input.Get("residual_vertices_y");
  TH1F * _residual_vertices_z = (TH1F*)input.Get("residual_vertices_z");
  TH1F * _residual_vertices_phi = (TH1F*)input.Get("residual_vertices_phi");
  TH1F * _residual_vertices_s = (TH1F*)input.Get("residual_vertices_s");
  TH1F * _residual_vertices_t = (TH1F*)input.Get("residual_vertices_t");

  _residual_vertices_x->GetXaxis()->SetTitle("x [cm]");
  _residual_vertices_x->GetYaxis()->SetTitle("events");
  _residual_vertices_y->GetXaxis()->SetTitle("y [cm]");
  _residual_vertices_y->GetYaxis()->SetTitle("events");
  _residual_vertices_z->GetXaxis()->SetTitle("z [cm]");
  _residual_vertices_z->GetYaxis()->SetTitle("events");
  _residual_vertices_t->GetYaxis()->SetTitle("events");

  TCanvas * c = new TCanvas("c","c",800,600);
  gStyle->SetOptStat(0);

  _vertices_x->Draw("e");
  c->Update();
  c->Print("vertices_x.png");

  _vertices_y->Draw("e");
  c->Update();
  c->Print("vertices_y.png");

  _vertices_phi->Draw("e");
  c->Update();
  c->Print("vertices_phi.png");

  _vertices_z->Draw("e");
  c->Update();
  c->Print("vertices_z.png");

  _vertices_s->Draw("e");
  c->Update();
  c->Print("vertices_s.png");

  _vertices_y_vs_x->Draw("colz");
  c->Update();
  c->Print("vertices_y_vs_x.png");

  _vertices_z_vs_phi->Draw("colz");
  c->Update();
  c->Print("vertices_z_vs_phi.png");


  _true_vertices_x->Draw("e");
  c->Update();
  c->Print("true_vertices_x.png");

  _true_vertices_y->Draw("e");
  c->Update();
  c->Print("true_vertices_y.png");

  _true_vertices_phi->Draw("e");
  c->Update();
  c->Print("true_vertices_phi.png");

  _true_vertices_z->Draw("e");
  c->Update();
  c->Print("true_vertices_z.png");

  _true_vertices_s->Draw("e");
  c->Update();
  c->Print("true_vertices_s.png");

  _true_vertices_y_vs_x->Draw("colz");
  c->Update();
  c->Print("true_vertices_y_vs_x.png");

  _true_vertices_z_vs_phi->Draw("colz");
  c->Update();
  c->Print("true_vertices_z_vs_phi.png");


  _ratio_vertices_x->Draw("e");
  c->Update();
  c->Print("ratio_vertices_x.png");

  _ratio_vertices_y->Draw("e");
  c->Update();
  c->Print("ratio_vertices_y.png");

  _ratio_vertices_phi->Draw("e");
  c->Update();
  c->Print("ratio_vertices_phi.png");

  _ratio_vertices_z->Draw("e");
  c->Update();
  c->Print("ratio_vertices_z.png");

  _ratio_vertices_s->Draw("e");
  c->Update();
  c->Print("ratio_vertices_s.png");

  _ratio_vertices_y_vs_x->Draw("colz");
  c->Update();
  c->Print("ratio_vertices_y_vs_x.png");

  _ratio_vertices_z_vs_phi->Draw("colz");
  c->Update();
  c->Print("ratio_vertices_z_vs_phi.png");

  _residual_vertices_s_t->Draw("colz");
  c->Update();
  c->Print("residual_vertices_s_t.png");

  _residual_vertices_s__vs__s->Draw("colz");
  c->Update();
  c->Print("residual_vertices_s__vs__s.png");



  do_fit(distance, _residual_vertices_x);
  c->Update();
  c->Print("residual_vertices_x.png");

  do_fit(distance, _residual_vertices_y);
  c->Update();
  c->Print("residual_vertices_y.png");

  do_fit(distance, _residual_vertices_z);
  c->Update();
  c->Print("residual_vertices_z.png");

  _residual_vertices_t->Draw("e");
  c->Update();
  c->Print("residual_vertices_t.png");

  do_fit(distance, _residual_vertices_s);
  c->Update();
  c->Print("residual_vertices_s.png");

  _residual_vertices_phi->Draw("e");
  c->Update();
  c->Print("residual_vertices_phi.png");


  input.Close();

  //exit(0);
  return 1;
}

void do_fit(double distance, TH1F * h){

  TFitResultPtr r;
  double A, mean, sigma, smean, ssigma;
  double A2, mean2, sigma2, smean2, ssigma2;

  double nsigma_one_gaus = 2.;
  double nsigma_two_gaus = 20.;

  double duration = 2e-6 * 1.e4; // (s)

  TF1 *one_gaussian = new TF1("one_gaussian", "gaus", -nsigma_one_gaus*distance,nsigma_one_gaus*distance);
  one_gaussian->SetParameters(1.,0.,distance/2.);
  one_gaussian->SetParLimits(2,0.,nsigma_two_gaus*distance);
  r = h->Fit("one_gaussian","S","",-nsigma_one_gaus*distance,nsigma_one_gaus*distance);
  A = r->Parameter(0);
  mean = r->Parameter(1);
  sigma = r->Parameter(2);

  TF1 *two_gaussians = new TF1("two_gaussians", "gaus(0) + gaus(3)", -nsigma_two_gaus*distance,nsigma_two_gaus*distance);

  two_gaussians->SetParameters(A,mean,sigma,1.,0.,nsigma_two_gaus*distance);
  two_gaussians->SetParLimits(0,0.5*A,2*A);
  two_gaussians->SetParLimits(1,0.5*mean, 2.*mean);
  two_gaussians->SetParLimits(2,0.5*sigma, 2.*sigma);
  two_gaussians->SetParLimits(5, nsigma_one_gaus*sigma, nsigma_two_gaus*distance);
  r = h->Fit("two_gaussians","S","",-nsigma_two_gaus*distance,nsigma_two_gaus*distance);
  A = r->Parameter(0);
  mean = r->Parameter(1);
  sigma = r->Parameter(2);
  smean = r->ParError(1);
  ssigma = r->ParError(2);
  A2 = r->Parameter(3);
  mean2 = r->Parameter(4);
  sigma2 = r->Parameter(5);
  smean2 = r->ParError(4);
  ssigma2 = r->ParError(5);

  TF1 *first_gaussian = new TF1("first_gaussian", "gaus", -nsigma_two_gaus*distance,nsigma_two_gaus*distance);
  first_gaussian->SetParameters(A, mean, sigma);

  TF1 *second_gaussian = new TF1("second_gaussian", "gaus", -nsigma_two_gaus*distance,nsigma_two_gaus*distance);
  second_gaussian->SetParameters(A2, mean2, sigma2);

  double I1 = first_gaussian->Integral(-nsigma_two_gaus*distance,nsigma_two_gaus*distance)/h->GetBinWidth(1);
  double I2 = second_gaussian->Integral(-nsigma_two_gaus*distance,nsigma_two_gaus*distance)/h->GetBinWidth(1);

  clog << " I1 " << I1 << " I2 " << I2 << " I1 + I2 " << I1 + I2 << endl;

  double rate2 = I2/duration;

  TPaveText *p_residual_vertices_x = new TPaveText(.6,.4,.9,.9,"NDC");
  p_residual_vertices_x->SetFillColor(0);
  p_residual_vertices_x->AddText(Form("mean = %.0f #pm %.0f", mean, smean));
  p_residual_vertices_x->AddText(Form("sigma = %.0f #pm %.0f", sigma, ssigma));
  p_residual_vertices_x->AddText(Form("sigma2 = %.0f #pm %.0f", sigma2, ssigma2));
  //  p_residual_vertices_x->AddText(Form("I1/I2 = %.2f", I1/I2));
  p_residual_vertices_x->AddText(Form("bkg rate = %.0f kHz", rate2/1.e3));
  h->Draw("e");
  p_residual_vertices_x->Draw("same");

  h->GetXaxis()->SetRangeUser(-2.5*sigma2, 2.5*sigma2);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(0.6);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(0.7);
  h->SetTitle("");

  return;
}
