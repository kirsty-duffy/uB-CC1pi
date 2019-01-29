#include "CC1pi_treevars.h"
#include "TMVA/Reader.h"

struct mupiBDT{
   float length_over_startend; // Track length vs start-end distance (division)
   float track_length; // Track length
   float ndaughters; // No. reco daughters
   float lnLmipovermu; // MIP vs mu LLH (is there a Bragg peak?)
   float lnLmipoverpi; // MIP vs pi LLH (is there a Bragg peak? Is this duplicating above?)
   //float pandoraclassedastrack; // Pandora track/shower classification (watch out for problems with DIC and the fact that this is really a bool)
   float perc_used_hits; // Percentage of used hits in track
   float residuals_mean; // Mean of track residuals (watch out for problems with DIC)
   float residuals_stddev; // Std dev of track residuals (watch out for problems with DIC)
   float MCS_pi_maxScatter; // Largest MCS scattering angle
   float MCS_pi_meanScatter; // Mean MCS scattering angle
   // TODO Energy near end of track not associated with a track

   std::string BookMVAType = "BDTG";
   std::string BookMVALoc = "/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2019-01-02/BDTtesting/dataset_vtxtrackprecut/weights/TMVAClassification_BDTG.weights.xml";
   TMVA::Reader *fReader;

   bool CheckMVAvars();
   void initialise_BDT();
};
