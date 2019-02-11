#include "CC1pi_treevars.h"
#include "TMVA/Reader.h"

struct mupiBDT{
   float length_over_startend; // Track length vs start-end distance (division)
   float track_length; // Track length
   float track_length_over_longest; // Track length divided by length of longest MIP track in event
   float ndaughters; // No. reco daughters
   float lnLmipovermu; // MIP vs mu LLH (is there a Bragg peak?)
   float lnLmipoverpi; // MIP vs pi LLH (is there a Bragg peak? Is this duplicating above?)
   float pandoraclassedastrack; // Pandora track/shower classification (watch out for problems with DIC and the fact that this is really a bool)
   float perc_used_hits; // Percentage of used hits in track
   float residuals_mean; // Mean of track residuals (watch out for problems with DIC)
   float residuals_stddev; // Std dev of track residuals (watch out for problems with DIC)
   float MCS_pi_maxScatter; // Largest MCS scattering angle
   float MCS_pi_meanScatter; // Mean MCS scattering angle
   float n_unused_hits_nearend; // No. hits near the end of the track not matched to any other track
   float unmatched_charge_nearend_plane2; // Total charge of hits on plane 2 near end of track not matched to any other track

   std::string BookMVAType_contained = "BDTG";
   std::string BookMVALoc_contained = "/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2019-01-28/dataset_cont/weights/TMVAClassification_cont_BDTG.weights.xml";
   TMVA::Reader *fReader_contained;
   std::string BookMVAType_exiting = "BDTG";
   std::string BookMVALoc_exiting = "/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2019-01-28/dataset_exiting/weights/TMVAClassification_exiting_BDTG.weights.xml";
   TMVA::Reader *fReader_exiting;

   bool CheckMVAvars();
   void initialise_BDT_contained();
   void initialise_BDT_exiting();
};
