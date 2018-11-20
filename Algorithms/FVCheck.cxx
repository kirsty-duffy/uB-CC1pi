#ifndef FVCHECK_CXX
#define FVCHECK_CXX

#include "FVCheck.h"

//TPC boundary + FV cut
const double FVxmin = 0 + 12;
const double FVxmax = 256.35 - 12;
const double FVymin = -115.53 + 35;
const double FVymax = 117.47 - 35;
const double FVzmin = 0.1 + 25;
const double FVzmax = 1036.9 - 85;

//Cut out additional section of mostly dead wires in z
const double deadzmin = 675.1;
const double deadzmax = 775.1;

bool inFV(double x, double y, double z){
   if (x > FVxmin && x < FVxmax && y > FVymin && y < FVymax && z > FVzmin && z < FVzmax && (z < deadzmin || z > deadzmax)) return true;
   else return false;
}

bool FVCheck(art::Event &evt, InputTags *CC1piInputTags)
{
   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel(CC1piInputTags->fSelectionLabel, selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "[FVCheck] SelectionResult product not found." << std::endl;
      throw std::exception();
   }
   std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
   art::fill_ptr_vector(selection_v, selection_h);

   // Get TPCObject
   art::FindManyP<ubana::TPCObject> tpcobject_from_selection(selection_h, evt, CC1piInputTags->fSelectionLabel);
   art::Ptr<ubana::TPCObject> tpcobj_candidate = tpcobject_from_selection.at(0).at(0);

   // Get neutrino candidate vertex from TPCObject
   recob::Vertex tpcobj_nu_vtx = tpcobj_candidate->GetVertex();
   double reco_nu_vtx[3];
   tpcobj_nu_vtx.XYZ(reco_nu_vtx);

   // space charge correction
   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
   std::vector<double> sce_corr = SCE->GetPosOffsets(reco_nu_vtx[0], reco_nu_vtx[1], reco_nu_vtx[2]);

   return inFV(reco_nu_vtx[0] + sce_corr.at(0), reco_nu_vtx[1] - sce_corr.at(1), reco_nu_vtx[2] - sce_corr.at(2));
}

#endif
