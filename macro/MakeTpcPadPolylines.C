// MakeTpcPadPolylines.C
//
// Standalone ROOT/sPHENIX macro that creates a ROOT file containing
// TPolyLine3D objects for TPC pad/timebin lines.
//
// Each saved polyline corresponds to one fixed (side, sector, layer, pad-in-sector).
// Its points are (x,y,z) for all requested time bins.
// x,y come from the TPC_FEE_CHANNEL_MAP-derived geometry convention.
// z is calculated from tbin using vdrift and t0.
//
// Run inside an sPHENIX environment, for example:
//   root -l -b -q 'MakeTpcPadPolylines.C(75570,"tpc_pad_polylines.root")'
//
// Useful options:
//   root -l -b -q 'MakeTpcPadPolylines.C(75570,"out.root",0.0076101,0.0,-1,-1,true,1,8)'
//     draw_pad_stride=1 draws all pads on canvases; this can be very heavy.
//     save_all_polylines=true saves every polyline object to the ROOT file.

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TPolyLine3D.h>
#include <TString.h>
#include <TSystem.h>
#include <TView.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <TParameter.h>
#include <phool/recoConsts.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace TpcPadPolylineInternal
{
  static const int N_SIDES = 2;
  static const int N_SECTORS = 12;
  static const int N_REGIONS = 3;
  static const int N_ROWS = 16;
  static const int N_LAYERS = 48;
  static const int FIRST_LAYER = 7;
  static const int LAST_LAYER = 54;
  static const int N_FEE = 26;
  static const int N_CH = 256;

  double wrap_phi(double phi)
  {
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    while (phi >   M_PI) phi -= 2.0 * M_PI;
    return phi;
  }

  int layer_to_region(const int layer)
  {
    if (layer < FIRST_LAYER || layer > LAST_LAYER) return -1;
    return (layer - FIRST_LAYER) / N_ROWS;
  }

  int pads_per_sector(const int region)
  {
    if (region == 0) return 96;
    if (region == 1) return 128;
    if (region == 2) return 192;
    return 0;
  }

  struct PadGeomFromCdb
  {
    double layer_radius_cm[N_LAYERS];
    double phi_bin_width_layer[N_LAYERS];
    double sec_min_phi[N_SIDES][N_SECTORS][N_REGIONS];
    double sec_max_phi[N_SIDES][N_SECTORS][N_REGIONS];

    PadGeomFromCdb()
    {
      for (int il = 0; il < N_LAYERS; ++il)
      {
        layer_radius_cm[il] = 0.0;
        phi_bin_width_layer[il] = 0.0;
      }
      for (int is = 0; is < N_SIDES; ++is)
      {
        for (int sec = 0; sec < N_SECTORS; ++sec)
        {
          for (int reg = 0; reg < N_REGIONS; ++reg)
          {
            sec_min_phi[is][sec][reg] = 0.0;
            sec_max_phi[is][sec][reg] = 0.0;
          }
        }
      }
    }

    void initialize_sector_phi_from_layer(const int layer_index,
                                          const double min_phi,
                                          const double max_phi,
                                          const double phi_bin_width)
    {
      if (layer_index < 0 || layer_index >= N_LAYERS) return;
      const int region = layer_index / N_ROWS;
      if (region < 0 || region >= N_REGIONS) return;

      // Same as PHG4/TpcPadMap convention.
      const double sector_phi = std::fabs(max_phi - min_phi) + phi_bin_width;

      for (int side = 0; side < N_SIDES; ++side)
      {
        for (int sector = 0; sector < N_SECTORS; ++sector)
        {
          const double sector_anchor = M_PI - 2.0 * M_PI / 12.0 * static_cast<double>(sector + 1);

          if (side == 0)
          {
            sec_min_phi[side][sector][region] = sector_anchor + (-max_phi - phi_bin_width / 2.0);
            sec_max_phi[side][sector][region] = sec_min_phi[side][sector][region] + sector_phi;
          }
          else
          {
            sec_max_phi[side][sector][region] = sector_anchor + (max_phi + phi_bin_width / 2.0);
            sec_min_phi[side][sector][region] = sec_max_phi[side][sector][region] - sector_phi;
          }
        }
      }
    }

    bool load_from_cdb(const int runnumber, const int verbosity)
    {
      // This is enough for CDBInterface::getUrl("TPC_FEE_CHANNEL_MAP") to pick run-dependent CDB.
      if (runnumber > 0)
      {
        recoConsts::instance()->set_IntFlag("RUNNUMBER", runnumber);
        CDBInterface::instance()->Verbosity(verbosity);
        recoConsts::instance()->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
        recoConsts::instance()->set_uint64Flag("TIMESTAMP", runnumber);
      }

      CDBInterface* cdb = CDBInterface::instance();
      const std::string calibdir = cdb->getUrl("TPC_FEE_CHANNEL_MAP");

      if (calibdir.empty())
      {
        std::cout << "ERROR: no TPC_FEE_CHANNEL_MAP found from CDBInterface" << std::endl;
        return false;
      }

      if (verbosity > 0)
      {
        std::cout << "Using TPC_FEE_CHANNEL_MAP: " << calibdir << std::endl;
      }

      CDBTTree cdbttree(calibdir);
      cdbttree.LoadCalibrations();

      std::vector<double> pad_phi[N_LAYERS];
      std::vector<double> pad_R[N_LAYERS];

      for (int fee = 0; fee < N_FEE; ++fee)
      {
        for (int ch = 0; ch < N_CH; ++ch)
        {
          const unsigned int key = 256U * static_cast<unsigned int>(fee) + static_cast<unsigned int>(ch);

          const int layer = cdbttree.GetIntValue(key, "layer");
          if (layer > 6 && layer < 55)
          {
            const int layer_index = layer - FIRST_LAYER;
            pad_phi[layer_index].push_back(cdbttree.GetDoubleValue(key, "phi"));
            pad_R[layer_index].push_back(cdbttree.GetDoubleValue(key, "R"));
          }
        }
      }

      for (int layer_index = 0; layer_index < N_LAYERS; ++layer_index)
      {
        if (pad_phi[layer_index].empty() || pad_R[layer_index].empty())
        {
          std::cout << "ERROR: missing pad map entries for layer "
                    << layer_index + FIRST_LAYER << std::endl;
          return false;
        }

        double radius_cm = 0.0;
        for (std::size_t i = 0; i < pad_R[layer_index].size(); ++i)
        {
          radius_cm += pad_R[layer_index][i];
        }
        radius_cm /= static_cast<double>(pad_R[layer_index].size());
        radius_cm /= 10.0;  // CDB R is mm; use cm.
        layer_radius_cm[layer_index] = radius_cm;

        double min_phi = *std::min_element(pad_phi[layer_index].begin(), pad_phi[layer_index].end());
        double max_phi = *std::max_element(pad_phi[layer_index].begin(), pad_phi[layer_index].end());

        // Match PHG4TpcDetector convention used in your TpcPadMap helper.
        min_phi -= M_PI / 2.0;
        max_phi -= M_PI / 2.0;

        const int layer = layer_index + FIRST_LAYER;
        const int region = layer_to_region(layer);
        const int n_pad_sec = pads_per_sector(region);

        if (n_pad_sec < 2)
        {
          std::cout << "ERROR: bad pads_per_sector for region " << region << std::endl;
          return false;
        }

        const double phi_bin_width = std::fabs(max_phi - min_phi) /
                                     static_cast<double>(n_pad_sec - 1);
        phi_bin_width_layer[layer_index] = phi_bin_width;

        initialize_sector_phi_from_layer(layer_index, min_phi, max_phi, phi_bin_width);
      }

      return true;
    }

    double global_phi(const int side,
                      const int sector,
                      const int layer,
                      const int pad_in_sector) const
    {
      const int region = layer_to_region(layer);
      const int layer_index = layer - FIRST_LAYER;
      const double phistep = phi_bin_width_layer[layer_index];

      // Center of pad, matching PHG4TpcCylinderGeom-style sector_max - (pad+0.5)*step.
      double phi = 0.0;
      if (side == 0)
      {
        phi = sec_min_phi[side][sector][region] +
              (static_cast<double>(pad_in_sector) + 0.5) * phistep;
      }
      else
      {
        phi = sec_max_phi[side][sector][region] -
              (static_cast<double>(pad_in_sector) + 0.5) * phistep;
      }
      return wrap_phi(phi);
    }
  };

  int default_ntbins_from_run(const int runnumber,
                             const double vdrift_cm_per_ns,
                             const double maxdriftlength_cm,
                             const int user_ntbins)
  {
    if (user_ntbins > 0) return user_ntbins;

    // Defaults copied from the logic in the user message.
    // Exact RunnumberRange constants are not available in this standalone macro,
    // so this uses practical switches and lets you override NTBins directly.
    double tpc_adc_clock_ns = 53.326184;
    double extended_readout_time_ns = 0.0;

    // Most Run-2 pp / early Run-3 20 MHz TPC clock.
    // Override tpc_adc_clock_ns explicitly if your run is in a different period.
    if (runnumber > 0)
    {
      tpc_adc_clock_ns = 53.326184;
    }

    const double minT = 0.0;
    const double maxT = extended_readout_time_ns + 2.0 * maxdriftlength_cm / vdrift_cm_per_ns;
    return static_cast<int>((maxT - minT) / tpc_adc_clock_ns) + 1;
  }

  double default_clock_from_mode(const int clock_mode)
  {
    // clock_mode: 0 = simulation / old default, 1 = 20 MHz, 2 = 17.5 MHz
    if (clock_mode == 1) return 50.037280;
    if (clock_mode == 2) return 56.881262;
    return 53.326184;
  }

  double default_extended_readout_from_mode(const int readout_mode)
  {
    // readout_mode: 0 = no extended readout, 1 = 24800 ns extended readout.
    if (readout_mode == 1) return 24800.0;
    return 0.0;
  }

  double z_from_timebin(const int side,
                        const int tbin,
                        const double tpc_adc_clock_ns,
                        const double vdrift_cm_per_ns,
                        const double t0_ns,
                        const double maxdriftlength_cm)
  {
    // Time coordinate for bin center.
    const double t_ns = (static_cast<double>(tbin) + 0.5) * tpc_adc_clock_ns;

    // Drift distance from the central membrane toward the endcap.
    double drift_cm = (t_ns - t0_ns) * vdrift_cm_per_ns;
    if (drift_cm < 0.0) drift_cm = 0.0;
    if (drift_cm > maxdriftlength_cm) drift_cm = maxdriftlength_cm;

    // Convention: side 0 negative z, side 1 positive z.
    // If your local convention is opposite, swap the signs here.
    return (side == 0) ? -drift_cm : drift_cm;
  }
}

void MakeTpcPadPolylines(const int runnumber = 75570,
                         const char* outfilename = "tpc_pad_polylines.root",
                         const double vdrift_cm_per_ns = 0.0076101,
                         const double t0_ns = 0.0,
                         const int ntbins_user = -1,
                         const int clock_mode = 0,
                         const bool save_all_polylines = true,
                         const int draw_pad_stride = 16,
                         const int draw_tbin_stride = 8,
                         const double maxdriftlength_cm = 102.325,
                         const int readout_mode = 0,
                         const int verbosity = 1)
{
  using namespace TpcPadPolylineInternal;

  // Needed for recoConsts if not already loaded by rootlogon.
  gSystem->Load("libffamodules.so");
  gSystem->Load("libcdbobjects.so");

  PadGeomFromCdb geom;
  if (!geom.load_from_cdb(runnumber, verbosity))
  {
    std::cout << "Failed to load CDB pad map. Stop." << std::endl;
    return;
  }

  const double tpc_adc_clock_ns = default_clock_from_mode(clock_mode);
  const double extended_readout_time_ns = default_extended_readout_from_mode(readout_mode);

  int ntbins = ntbins_user;
  if (ntbins <= 0)
  {
    const double maxT = extended_readout_time_ns + 2.0 * maxdriftlength_cm / vdrift_cm_per_ns;
    ntbins = static_cast<int>(maxT / tpc_adc_clock_ns) + 1;
  }


  std::cout << "Output file: " << outfilename << std::endl;
  std::cout << "vdrift = " << vdrift_cm_per_ns << " cm/ns" << std::endl;
  std::cout << "t0     = " << t0_ns << " ns" << std::endl;
  std::cout << "clock  = " << tpc_adc_clock_ns << " ns" << std::endl;
  std::cout << "extended_readout_time = " << extended_readout_time_ns << " ns" << std::endl;
  std::cout << "maxdriftlength = " << maxdriftlength_cm << " cm" << std::endl;
  std::cout << "NTBins = " << ntbins << std::endl;

  TFile* fout = new TFile(outfilename, "RECREATE");
  fout->cd();

  TDirectory* dmeta = fout->mkdir("metadata");
  dmeta->cd();
  TNamed meta_run("runnumber", Form("%d", runnumber));
  TNamed meta_vd("vdrift_cm_per_ns", Form("%.12g", vdrift_cm_per_ns));
  TNamed meta_t0("t0_ns", Form("%.12g", t0_ns));
  TNamed meta_clock("tpc_adc_clock_ns", Form("%.12g", tpc_adc_clock_ns));
  TNamed meta_ntbins("NTBins", Form("%d", ntbins));
  meta_run.Write();
  meta_vd.Write();
  meta_t0.Write();
  meta_clock.Write();
  meta_ntbins.Write();

  TDirectory* dpoly = fout->mkdir("polylines");
  TDirectory* dcanv = fout->mkdir("canvases");

  Long64_t n_saved = 0;
  Long64_t n_drawn = 0;

  for (int side = 0; side < N_SIDES; ++side)
  {
    for (int sector = 0; sector < N_SECTORS; ++sector)
    {
      for (int region = 0; region < N_REGIONS; ++region)
      {
        const int first_layer = FIRST_LAYER + region * N_ROWS;
        const int last_layer = first_layer + N_ROWS - 1;
        const int n_pads = pads_per_sector(region);

        TString cname = Form("c_pad_polylines_side%d_sector%02d_region%d", side, sector, region);
        TCanvas* c = new TCanvas(cname, cname, 1200, 900);
        c->cd();
        TView* view = TView::CreateView(1);
        view->SetRange(-85.0, -85.0, -maxdriftlength_cm,
                        85.0,  85.0,  maxdriftlength_cm);

        for (int layer = first_layer; layer <= last_layer; ++layer)
        {
          const int layer_index = layer - FIRST_LAYER;
          const double r_cm = geom.layer_radius_cm[layer_index];

          for (int pad = 0; pad < n_pads; ++pad)
          {
            const double phi = geom.global_phi(side, sector, layer, pad);
            const double x_cm = r_cm * std::cos(phi);
            const double y_cm = r_cm * std::sin(phi);

            TPolyLine3D* pl = new TPolyLine3D(ntbins);
            //pl->SetName(Form("pl_s%d_sec%02d_reg%d_l%02d_pad%03d", side, sector, region, layer, pad));
            //pl->SetTitle(Form("side=%d sector=%d region=%d layer=%d pad=%d;X [cm];Y [cm];Z [cm]", side, sector, region, layer, pad));
            pl->SetLineWidth(1);

            for (int it = 0; it < ntbins; ++it)
            {
              const double z_cm = z_from_timebin(side, it, tpc_adc_clock_ns,
                                                 vdrift_cm_per_ns, t0_ns,
                                                 maxdriftlength_cm);
              pl->SetPoint(it, x_cm, y_cm, z_cm);
            }

            if (save_all_polylines)
            {
              dpoly->cd();
              const TString plname = Form("pl_s%d_sec%02d_reg%d_l%02d_pad%03d", side, sector, region, layer, pad);
              pl->Write(plname);
              ++n_saved;
            }

            // Draw a controlled subset by default. Set draw_pad_stride=1 to draw every pad.
            if (draw_pad_stride > 0 && draw_tbin_stride > 0 && (pad % draw_pad_stride) == 0)
            {
              // Draw a thinned visual copy if requested. The saved polyline above remains full NTBins.
              const int nvis = ntbins / draw_tbin_stride + 1;
              TPolyLine3D* plvis = new TPolyLine3D(nvis);
              //plvis->SetName(Form("vis_%s", pl->GetName()));
              int iv = 0;
              for (int it = 0; it < ntbins; it += draw_tbin_stride)
              {
                const double z_cm = z_from_timebin(side, it, tpc_adc_clock_ns,
                                                   vdrift_cm_per_ns, t0_ns,
                                                   maxdriftlength_cm);
                plvis->SetPoint(iv, x_cm, y_cm, z_cm);
                ++iv;
              }
              plvis->Draw((n_drawn == 0) ? "" : "same");
              ++n_drawn;
            }

            if (!save_all_polylines)
            {
              delete pl;
            }
          }
        }

        dcanv->cd();
        c->Write();
        delete c;
      }
    }
  }

  fout->cd();
  fout->Write();
  fout->Close();

  std::cout << "Done." << std::endl;
  std::cout << "Saved full TPolyLine3D objects: " << n_saved << std::endl;
  std::cout << "Drawn visual polylines on canvases: " << n_drawn << std::endl;
}
