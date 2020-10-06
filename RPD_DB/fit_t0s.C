/*
 * fit_t0s.C
 *
 *  Created on: Feb 8, 2010
 *      Author: Promme
 *
 *      fitting t0 calibration constants by using the histograms
 *      created by RPD::DecodeChipDigits() svn vers. 190
 *
 *      If you use preselected mDST containing several runs, you might
 *      use RPD::Calibrate() instead for a faster solution
 */

#include "TH1F.h"
#include <sstream>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TBrowser.h"
#include <fstream>

void fit_t0s(){
	gROOT->SetStyle("Plain");
	const int NSLATS[2]={12,24};
	TFile* file = new TFile("hist.root");
	if (file->IsZombie()) return;
	file->cd();
	TH1F* hist_t0_calib[2][24][2];
	for (int iring = 0; iring < 2; iring++){
		for (int islat = 0; islat < NSLATS[iring]; islat++){
			for (int ipm = 0; ipm < 2; ipm++){
				stringstream title;
				title << "Vertex_in_targetregion/hist_t0_calib_" << iring << "_" << islat << "_" << ipm << endl;
				hist_t0_calib[iring][islat][ipm] = (TH1F*)file->Get(title.str().c_str());
			}
		}
	}
	double t0_calib_const_mean[2][24][2];
	double t0_calib_const_sigma[2][24][2];
	// draw, fit and write out the t0 calibration values
	TCanvas* t0_window = new TCanvas("t0_window", "t0s", 1200, 800);
	TCanvas* t0_window_one  = new TCanvas("t0_window_one", "t0", 800, 300);
	t0_window->Divide(12,6);
	int ipad(0);
	ofstream t0_calib_file;
	t0_calib_file.open ("fitted_t0_calib_constants.txt");

	for (int iring = 0; iring < 2; iring++){
		for (int islat = 0; islat < NSLATS[iring]; islat++){
			for (int ipm = 0; ipm < 2; ipm++){
				ipad++;
				t0_window->cd(ipad);
				gPad->Clear();
				if (hist_t0_calib[iring][islat][ipm]){
					hist_t0_calib[iring][islat][ipm]->Draw();
					stringstream title;
					title << "gaus_" << iring << "_" << islat << "_" << ipm << endl;
					float mean = hist_t0_calib[iring][islat][ipm]->GetMean();
					float maximum = hist_t0_calib[iring][islat][ipm]->GetBinCenter(hist_t0_calib[iring][islat][ipm]->GetMaximumBin());
					float rms  = hist_t0_calib[iring][islat][ipm]->GetRMS();
					TF1* gaus = new TF1(title.str().c_str(), "gaus(0)",maximum-5, maximum+5);
					gaus->SetParameters(hist_t0_calib[iring][islat][ipm]->GetMaximum()/rms,
							maximum,
							rms);
					hist_t0_calib[iring][islat][ipm]->Fit(gaus, "RB");
					gaus = hist_t0_calib[iring][islat][ipm]->GetFunction(title.str().c_str());
					t0_calib_file<< - gaus->GetParameter(1) << ", ";
					t0_calib_const_mean[iring][islat][ipm] = gaus->GetParameter(1);
					t0_calib_const_sigma[iring][islat][ipm]= gaus->GetParameter(2);
					gPad->Update();
					t0_window_one->cd(0);
					gPad->Clear();
					t0_window_one->SetTitle(title.str().c_str());
					hist_t0_calib[iring][islat][ipm]->Draw();
					gPad->Update();
					//sleep(1);
				}
			}
		}
		t0_calib_file << endl;
	}
	// calculate the mean value for ring A|B pm up|down
	double t0_mean[2][2]={{0.,0.},{0.,0.}};
	double t0_sum [2][2]={{0.,0.},{0.,0.}};
	for (int iring = 0; iring < 2; iring++){
		for (int islat = 0; islat < NSLATS[iring]; islat++){
			for (int ipm = 0; ipm < 2; ipm++){
				t0_mean[iring][ipm]+= t0_calib_const_mean[iring][islat][ipm];
				t0_sum[iring][ipm] += 1.;
			}
		}
	}
	for (int iring = 0; iring < 2; iring++){
		for (int ipm = 0; ipm < 2; ipm++){
			t0_mean[iring][ipm]/= t0_sum[iring][ipm];
			t0_calib_file << " mean value for ring " << iring << " pm " << ipm << " " << t0_mean[iring][ipm] << endl;
		}
	}
	// calculate the rms value
	double t0_rms [2][2]={{0.,0.},{0.,0.}};
	for (int iring = 0; iring < 2; iring++){
		for (int islat = 0; islat < NSLATS[iring]; islat++){
			for (int ipm = 0; ipm < 2; ipm++){
				t0_rms[iring][ipm]+= pow(t0_mean[iring][ipm]-t0_calib_const_mean[iring][islat][ipm],2);
			}
		}
	}
	for (int iring = 0; iring < 2; iring++){
		for (int ipm = 0; ipm < 2; ipm++){
			t0_rms[iring][ipm]/= pow(t0_sum[iring][ipm]+1,2);
			t0_rms[iring][ipm] = sqrt(t0_rms[iring][ipm]);
			t0_calib_file << " rms value for ring " << iring << " pm " << ipm << " " << t0_rms[iring][ipm] << endl;
		}
	}

	t0_calib_file.close();

	TBrowser *browser = new TBrowser();
}
/*
for (int iring = 0; iring < 2; iring++){
	for (int islat = 0; islat < NSLATS[iring]; islat++){
		for (int ipm = 0; ipm < 2; ipm++){
			//if ((iring == 0) && (islat >= 12)) continue;
			pos_const_t0++;
			Calib_constants_db[runnumb].t0[iring][islat][ipm] = const_t0[pos_const_t0];
			// test the filling procedure
			if (verbose){
				cout << " run number " << Calib_constants_db[runnumb].runnumb;
				cout << " ring " << iring << " slat " << islat;
				cout << " pm " << ipm;
				cout << " t0 const " << Calib_constants_db[runnumb].t0[iring][islat][ipm];
				cout << endl;
			}
		}
	}
}*/
