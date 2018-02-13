#include "FVCheck.h"

//TPC boundary + FV cut
const double FVxmin = 0 + 12;
const double FVxmax = 256.35 - 12;
const double FVymin = -115.53 + 35;
const double FVymax = 117.47 - 35;
const double FVzmin = 0.1 + 25;
const double FVzmax = 1036.9 - 85;

//Cut out additional section of mostly dead wires in z
const double deadzmin = 675;
const double deadzmax = 775;

bool inFV(double x, double y, double z){
   if (x > FVxmin && x < FVxmax && y > FVymin && y < FVymax && z > FVzmin && z < FVzmax && (z < deadzmin || z > deadzmax)) return true;
   else return false;
}

bool FVCheck(art::Event &evt)
{
   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel("UBXSec", selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found." << std::endl;
      throw std::exception();
   }
   std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
   art::fill_ptr_vector(selection_v, selection_h);

   // Get TPCObject
   art::FindManyP<ubana::TPCObject> tpcobject_from_selection(selection_h, evt, "UBXSec");
   art::Ptr<ubana::TPCObject> tpcobj_candidate = tpcobject_from_selection.at(0).at(0);

   // Get neutrino candidate vertex from TPCObject
   recob::Vertex tpcobj_nu_vtx = tpcobj_candidate->GetVertex();
   double reco_nu_vtx[3];
   tpcobj_nu_vtx.XYZ(reco_nu_vtx);

   return inFV(reco_nu_vtx[0], reco_nu_vtx[1], reco_nu_vtx[2]);
}
