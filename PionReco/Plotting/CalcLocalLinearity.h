#ifndef __CALCLOCALLINEARITY_H__
#define __CALCLOCALLINEARITY_H__
// Mean, covariance and std deviation functions for use below
// Copied from Marco: https://github.com/marcodeltutto/HitCosmicTag/blob/1bb389047776511774d3e9283985c8284879e1c9/Base/Tools.cxx
double mean(const std::vector<double>& data)
{
  if(data.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to mean" << std::endl;

  double result = 0.0;

  for(const auto& d : data)
    result += d;

  return (result / ((double)data.size()));
}

double cov (const std::vector<double>& data1,
            const std::vector<double>& data2,
            const std::vector<double>& data3)
{
  if(data1.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;
  if(data2.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;
  if(data3.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to cov" << std::endl;

  double result = 0.0;
  auto   mean1  = mean(data1);
  auto   mean2  = mean(data2);
  auto   mean3  = mean(data3);

  for(size_t i = 0; i < data1.size(); ++i)
    result += (data1[i] - mean1)*(data2[i] - mean2)*(data3[i] - mean3);

  return result/((double)data1.size());

}

double stdev(const std::vector<double>& data)
{
  if(data.size() == 0) std::cout << __PRETTY_FUNCTION__ << "You have me nill to stdev" << std::endl;

  double result = 0.0;
  auto    avg   = mean(data);
  for(const auto& d: data)
    result += (d - avg)*(d - avg);

  return std::sqrt(result/((double)data.size()));
}




// Function to calculate local linearity
// Copied from Marco: https://github.com/marcodeltutto/HitCosmicTag/blob/1bb389047776511774d3e9283985c8284879e1c9/Algorithms/ClassicLocalLinearityCalculator.cxx
/*
  std::vector<double> R;
  R.reserve(spacepoints.size());

  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  X.reserve(_slider_window);
  Y.reserve(_slider_window);
  Z.reserve(_slider_window);

  size_t n_either_side = std::ceil((_slider_window)/2);
  std::cout << "n_either_side = " << n_either_side << std::endl;

  for (size_t i_pos=0; i_pos<spacepoints.size(); i_pos++){
    // If there aren't enough hits before this hit in the vector, return 1 because we can't get a good measurement
    if (i_pos<n_either_side){
      R.push_back(0.0);
    }
    // If there aren't enough hits after this hit in the vector, return 1 because we can't get a good measurement
    else if ((spacepoints.size()-i_pos)<n_either_side){
      R.push_back(0.0);
    }
    else{ // this is the main function (for when we're looking at the middle of the track)
      // Fill X,Y,Z vectors
      std::cout << "Hit " << i_pos << ": Looking at " << i_pos-n_either_side << " to " << i_pos+n_either_side << std::endl;

      for (size_t i_w=i_pos-n_either_side; i_w<i_pos+n_either_side; i_w++){
        X.push_back(spacepoints.at(i_w).at(0));
        Y.push_back(spacepoints.at(i_w).at(1));
        Z.push_back(spacepoints.at(i_w).at(2));
      }

      auto c = cov(X,Y,Z);
      auto sX = stdev(X);
      auto sY = stdev(Y);
      auto sZ = stdev(Z);
      auto r = std::abs(c/(sX*sY*sZ));

      if (std::isnan(r)) r = 0.0;
      R.push_back(r);

      X.clear(); Y.clear(); Z.clear();
    }
  }
  return R;
};*/

// Another approach to "local linearity": calculate the dot product between the average vector from point to point before the point you're interested in and after that point. If linear, it should be about 1.
std::vector<double> GetLocalLinearityVec(std::vector<std::pair<std::vector<double>,int>> spacepoints, int _slider_window){

  std::vector<double> R;
  R.reserve(spacepoints.size());

  size_t n_either_side = std::ceil((_slider_window)/2);

  for (size_t i_pos=0; i_pos<spacepoints.size(); i_pos++){
    // If there aren't enough hits before this hit in the vector, return 1 because we can't get a good measurement
    if (i_pos<n_either_side){
      R.push_back(1.0);
    }
    // If there aren't enough hits after this hit in the vector, return 1 because we can't get a good measurement
    else if ((spacepoints.size()-i_pos-1)<n_either_side){
      R.push_back(1.0);
    }
    else{ // this is the main function (for when we're looking at the middle of the track)

      // Calculate average vector direction before the point we're interested in
      double theta_before=0;
      double phi_before=0;
      for (size_t i_w=i_pos-n_either_side; i_w<i_pos; i_w++){
        // First get the vector betwen this point and the next one
        TVector3 vec0(spacepoints.at(i_w).first.at(0),spacepoints.at(i_w).first.at(1),spacepoints.at(i_w).first.at(2));
        TVector3 vec1(spacepoints.at(i_w+1).first.at(0),spacepoints.at(i_w+1).first.at(1),spacepoints.at(i_w+1).first.at(2));

        TVector3 diff = vec1-vec0;
        theta_before += diff.Theta();
        phi_before += diff.Phi();
      }
      TVector3 avg_before(1,1,1);
      avg_before.SetMag(1);
      avg_before.SetTheta(theta_before/(double)n_either_side);
      avg_before.SetPhi(phi_before/(double)n_either_side);

      // Calculate average vector direction after the point we're interested in
      double theta_after=0;
      double phi_after=0;
      for (size_t i_w=i_pos; i_w<i_pos+n_either_side; i_w++){
        // First get the vector betwen this point and the next one
        TVector3 vec0(spacepoints.at(i_w).first.at(0),spacepoints.at(i_w).first.at(1),spacepoints.at(i_w).first.at(2));
        TVector3 vec1(spacepoints.at(i_w+1).first.at(0),spacepoints.at(i_w+1).first.at(1),spacepoints.at(i_w+1).first.at(2));

        TVector3 diff = vec1-vec0;
        theta_after += diff.Theta();
        phi_after += diff.Phi();
      }
      TVector3 avg_after(1,1,1);
      avg_after.SetMag(1);
      avg_after.SetTheta(theta_after/(double)n_either_side);
      avg_after.SetPhi(phi_after/(double)n_either_side);

      // Now calculate the difference in angles between those two vectors
      double r = avg_after.Dot(avg_before);
      R.push_back(r);
    }
  }
  return R;
};

// For each spacepoint, calculate angle between track start direction and the direction at this space point (i.e. the direction of the vector between this spacepoint and the next)
std::vector<double> GetSPDirVec(std::vector<std::pair<std::vector<double>,int>> spacepoints, double track_theta=-9999, double track_phi=-9999){

  std::vector<double> dir;
  dir.reserve(spacepoints.size());

  TVector3 beam_dir(1,1,1);
  if (track_theta==-9999 && track_phi==-9999){
    std::cout << "No track theta or phi given, setting to default beam direction" << std::endl;
    beam_dir = TVector3(0,0,1);
  }
  else{
    beam_dir.SetTheta(track_theta);
    beam_dir.SetPhi(track_phi);
  }
  beam_dir.SetMag(1.0);

  for (size_t i_pos=0; i_pos<spacepoints.size()-1; i_pos++){
      // For every point, store the direction vector between this point and the next
      // Won't work for the last point, so for that one just store the same as the previous one

      // First get the vector betwen this point and the next one
      TVector3 vec0(spacepoints.at(i_pos).first.at(0),spacepoints.at(i_pos).first.at(1),spacepoints.at(i_pos).first.at(2));
      TVector3 vec1(spacepoints.at(i_pos+1).first.at(0),spacepoints.at(i_pos+1).first.at(1),spacepoints.at(i_pos+1).first.at(2));

      TVector3 diff = vec1-vec0;
      diff.SetMag(1.0);

      // Calculate vector direction (direction w.r.t. the beam, which I assume is a vector (0,0,1)) for every point before the point we're interested in
      dir.push_back(TMath::ACos(beam_dir.Dot(diff)));
      }

      // Last point: repeat second-last point
      dir.push_back(dir.at(dir.size()-1));

  return dir;
};

// For each spacepoint, calculate angle between track start direction and the direction at this space point (i.e. the direction of the vector between this spacepoint and the next)
std::vector<double> GetCuSumVec(std::vector<double> dirs){

  std::vector<double> cs;
  cs.reserve(dirs.size());

  double sum = std::accumulate(std::begin(dirs), std::end(dirs), 0.0);
  double mean =  sum / dirs.size();

  double cusum = 0;
  for (size_t i_pos=0; i_pos<dirs.size(); i_pos++){
    cusum += (dirs.at(i_pos)-mean);
    cs.push_back(cusum);
  }

  return cs;
};

// Try calculating Welch's t-test, following https://www.nature.com/articles/nature03528.pdf and equation in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2134886/ and https://en.wikipedia.org/wiki/Welch%27s_t-test
std::vector<double> GetWelchttestVec(std::vector<double> dirs, int _slider_window=0){

  std::vector<double> R;
  R.reserve(dirs.size());

  size_t n_either_side = _slider_window+1;//std::ceil((_slider_window)/2);
  for (size_t i_pos=0; i_pos<dirs.size(); i_pos++){
    // If there aren't enough hits before this hit in the vector, return 1 because we can't get a good measurement
    if (i_pos<n_either_side){
      R.push_back(1.0);
    }
    // If there aren't enough hits after this hit in the vector, return 1 because we can't get a good measurement
    else if ((dirs.size()-i_pos-1)<n_either_side){
      R.push_back(1.0);
    }
    else{ // this is the main function (for when we're looking at the middle of the track)
      // Welch's t-test: t = (mean1-mean2)/(sqrt(var1/N1+var2/N2))

      int istart_before, iend_before, istart_after, iend_after;

      if (_slider_window==0){
        // If no sider window is given, split track at i_pos and use all hits before and all hits after
        istart_before = 0;
        iend_before = i_pos-1;

        istart_after = i_pos;
        iend_after = dirs.size()-1;
      }
      else{
        // If a slider window is given, split track and take n_either_side-1 hits on either side of i_pos
        istart_before = i_pos-n_either_side;
        iend_before = i_pos-1;

        istart_after = i_pos;
        iend_after = i_pos+n_either_side-1;
      }

      int n_before = iend_before-istart_before;
      int n_after = iend_after-istart_after;

      double sum_before = std::accumulate(std::begin(dirs)+istart_before, std::begin(dirs)+iend_before, 0.0);
      double mean_before =  sum_before / n_before;

      double accum = 0.0;
      std::for_each (std::begin(dirs)+istart_before, std::begin(dirs)+iend_before, [&](const double d) {
        accum += (d - mean_before) * (d - mean_before);
      });

      double var_before = accum / (n_before-1);

      double sum_after = std::accumulate(std::begin(dirs)+istart_after, std::begin(dirs)+iend_after, 0.0);
      double mean_after =  sum_after / n_after;

      accum = 0.0;
      std::for_each (std::begin(dirs)+istart_after, std::begin(dirs)+iend_after, [&](const double d) {
        accum += (d - mean_after) * (d - mean_after);
      });

      double var_after = accum / (n_after-1);

      double t = (mean_before-mean_after)/(TMath::Sqrt((var_before/n_before)+(var_after/n_after)));
      R.push_back(TMath::Abs(t));
    }
  }

  return R;
};



// Before we can do local linearity we have to make sure the spacepoints are correctly ordered
// A simple approach: start with the spacepoint closest to the track start (or reconstructed neutrino vertex?) and add the closest hit that hasn't been added yet until we have them all
std::vector<std::pair<std::vector<double>,int>> OrderSpacepoints(std::vector<std::vector<double>> spacepoints, std::vector<double> startpoint=std::vector<double>()){
  std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints;

  // Check input vector contains some points
  if (spacepoints.size()==0){
    std::cout << "Error: spacepoints vector is empty. Returning an empty vector." << std::endl;
    return ordered_spacepoints;
  }

  int firstpoint_idx = -9999;

  // If given a start point, find the space point closest to it
  if (startpoint.size()>0){
    //Check start point is valid
    if (startpoint.at(0)==-999 && startpoint.at(1)==-999 && startpoint.at(2)==-999){
      std::cout << "Error: cannot order spacepoints without a valid starting point. Returning an empty vector." << std::endl;
      return ordered_spacepoints;
    }

    // Find the closest spacepoint to the start point
    double dist = 9999;
    TVector3 start(startpoint.at(0),startpoint.at(1),startpoint.at(2));
    TVector3 firstpoint(1,1,1);
    for (int i=0; i<(int)spacepoints.size(); i++){
      std::vector<double> sp = spacepoints.at(i);
      TVector3 testpoint(sp.at(0),sp.at(1),sp.at(2));
      double testdist = (testpoint-start).Mag();
      if (testdist<dist){
        dist = testdist;
        firstpoint = testpoint;
        firstpoint_idx = i;
      }
    } // end loop to find closest spacepoint to the start
  }
  else{
    // Asking for the track start point doesn't seem to work well because Pandora won't necessarily give us what we want (which is a point at one end of the track - doesn't really matter which)
    // To save time, start out by finding the spacepoints with maximum and minimum x, y, and z. Those can be our candidates for the end points (may not always be exactly true, but probably a good start most of the time)
    std::vector<int> idx = {-999,-999,-999,-999,-999,-999};
    double xmin=99999, xmax=-99999, ymin=99999, ymax=-99999, zmin=99999, zmax=-99999;
    for(int i=0; i<(int)spacepoints.size(); i++){
      std::vector<double> sp = spacepoints.at(i);
      if (sp.at(0)<xmin){
        idx.at(0) = i;
        xmin = sp.at(0);
      }
      if (sp.at(0)>xmax){
        idx.at(1) = i;
        xmax = sp.at(0);
      }
      if (sp.at(1)<ymin){
        idx.at(2) = i;
        ymin = sp.at(1);
      }
      if (sp.at(1)>ymax){
        idx.at(3) = i;
        ymax = sp.at(1);
      }
      if (sp.at(2)<zmin){
        idx.at(4) = i;
        zmin = sp.at(2);
      }
      if (sp.at(2)>zmax){
        idx.at(5) = i;
        zmax = sp.at(2);
      }
    } // end loop over spacepoints to find six extremes

    // std::cout << xmin << "  " << idx.at(0) << std::endl;
    //   std::cout << xmax << "  " << idx.at(1) << std::endl;
    //     std::cout << ymin << "  " << idx.at(2) << std::endl;
    //       std::cout << ymax << "  " << idx.at(3) << std::endl;
    //         std::cout << zmin << "  " << idx.at(4) << std::endl;
    //           std::cout << zmax << "  " << idx.at(5) << std::endl;

    // Now examine these six extreme points. For each one, calculate the summed distance to all other spacepoints. The point with the maximum summed distance is our "most extreme" point = track end candidate
    std::vector<double> dist = {0,0,0,0,0,0};
    for (int i=0; i<(int)spacepoints.size(); i++){
      std::vector<double> sp = spacepoints.at(i);
      TVector3 testpoint(sp.at(0),sp.at(1),sp.at(2));

      for (int i_ext=0; i_ext<(int)idx.size(); i_ext++){
        if (idx.at(i_ext)==-999) continue;
        std::vector<double> extremepoint = spacepoints.at(idx.at(i_ext));
        TVector3 extpt(extremepoint.at(0), extremepoint.at(1), extremepoint.at(2));
        dist.at(i_ext) += ((extpt-testpoint).Mag());
      }
    } // end loop to find closest spacepoint to the start

    // Now find the maximum entry in dist, and use the corresponding entry in idx to find its index in spacepoints.
    firstpoint_idx = std::distance(dist.begin(), std::max_element(dist.begin(), dist.end()));
  } // end else (condition if no start point is given)

  // push back first point into new vector of spacepoints
  ordered_spacepoints.push_back(std::make_pair(spacepoints.at(firstpoint_idx),firstpoint_idx));
  std::vector<int> used_indices;
  used_indices.push_back(firstpoint_idx);

  // Now loop through all spacepoints. For each one, find the next closest spacepoint that hasn't been seen already and add it to ordered_spacepoints
  while (used_indices.size() < spacepoints.size()){

    std::vector<double> curr_sp = ordered_spacepoints.at(ordered_spacepoints.size()-1).first;
    TVector3 currentpoint(curr_sp.at(0),curr_sp.at(1),curr_sp.at(2));

    TVector3 nextpoint(1,1,1);
    double dist = 9999;
    int nextpoint_idx = -9999;
    for (int i=0; i<(int)spacepoints.size(); i++){
      // Check whether this spacepoint has already been used. If so, continue
      if (std::find(used_indices.begin(), used_indices.end(), i) != used_indices.end()) continue;

      // Now look for the closest spacepoint
      std::vector<double> sp = spacepoints.at(i);
      TVector3 testpoint(sp.at(0),sp.at(1),sp.at(2));
      double testdist = (testpoint-currentpoint).Mag();
      if (testdist<dist){
        dist = testdist;
        nextpoint = testpoint;
        nextpoint_idx = i;
      }
    } // end loop to find closest spacepoint to the current spacepoint

    // Now we have found the closest spacepoint, push it back in the vector
    std::vector<double> tmp = {nextpoint.X(),nextpoint.Y(),nextpoint.Z()};
    ordered_spacepoints.push_back(std::make_pair(tmp,nextpoint_idx));
    used_indices.push_back(nextpoint_idx);
  }

  return ordered_spacepoints;
};

int FindKinkIdx(std::vector<double> *dirs, std::vector<double> *cusum=nullptr, std::vector<double> *ttest=nullptr){
  // If cusum and ttest vectors are not given, calculate them!
  std::vector<double> cusum_tmp;
  std::vector<double> ttest_tmp;
  if (cusum==nullptr){
    cusum_tmp = GetCuSumVec(*dirs);
    cusum = &cusum_tmp;
  }
  if (ttest==nullptr){
    ttest_tmp = GetWelchttestVec(*dirs);
    ttest = &ttest_tmp;
  }

  // Loop through cusum vector and find point with maximum difference from 0
  std::vector<double>::iterator cusum_max_it = std::max_element(cusum->begin(),cusum->end());
  std::vector<double>::iterator cusum_min_it = std::min_element(cusum->begin(),cusum->end());

  int cusum_absmax_idx = std::distance(cusum->begin(),cusum_max_it);

  if (TMath::Abs(*cusum_min_it)>TMath::Abs(*cusum_max_it)) cusum_absmax_idx = std::distance(cusum->begin(),cusum_min_it);

  // Look at spacepoints near the one with maximum difference from 0 to find the one with maximum t-test value
  // "near" is tuneable! +/- 10 for now
  int n_near=10;
  int i_low = 0, i_high=ttest->size()-1;
  if (cusum_absmax_idx>=n_near) i_low = cusum_absmax_idx-n_near;
  if (cusum_absmax_idx+n_near<ttest->size()) i_high = cusum_absmax_idx+n_near;
  std::vector<double>::iterator ttest_max_it = std::max_element(ttest->begin()+i_low,ttest->begin()+i_high);

  // Check: is this maximum t-test value greater than threshold? If not, end the function as we have not found a decent kink
  // Threshold is tuneable! 5 for now
  if (*ttest_max_it<10){
    return -9999;
  }

  // If maximum t-test value is greater than threshold, we have a kink! Return its index so we can split vectors in two
  return std::distance(ttest->begin(),ttest_max_it);

};

std::vector<int>* SplitTracks(std::vector<double> dirs, std::vector<double> cusum = std::vector<double>(), std::vector<double> ttest = std::vector<double>()){
  // If cusum and ttest vectors are not given, calculate them!
  if (cusum.size()==0){
    cusum = GetCuSumVec(dirs);
  }
  if (ttest.size()==0){
    ttest = GetWelchttestVec(dirs);
  }

  // Make vectors of vectors so we can loop
  std::vector<std::vector<double>> dirs_vec;
  std::vector<std::vector<double>> cusum_vec;
  std::vector<std::vector<double>> ttest_vec;
  dirs_vec.push_back(dirs);
  cusum_vec.push_back(cusum);
  ttest_vec.push_back(ttest);

  // Make vector to store kink indices
  std::vector<int> *kink_indices = new std::vector<int>();

  int nrpts=0;
  while(nrpts<10){
    // this loop goes over the vectors of spacepoints we are considering.
    // It keeps going indefinitely until one of the break conditions is reached

    // Break conditions: no more points with max. cusum difference from 0 and t-test value above threshold, or maximum number of iterations has been reached, or track has been split to be too short
    // Maxumum number of iterations is tuneable! 10 for now
    // Minimum track length (in terms of no. spacepoints) is tuneable! 5 for now

    std::vector<std::vector<double>> dirs_vec_tmp;
    std::vector<std::vector<double>> cusum_vec_tmp;
    std::vector<std::vector<double>> ttest_vec_tmp;

    int kink_indices_size_before=kink_indices->size();

    for (int i=0; i<dirs_vec.size(); i++){
      // if (dirs_vec.at(i).size()<5) continue;

      int kinkidx = -9999;
      if (i==0 && nrpts==0){
        kinkidx = FindKinkIdx(&(dirs_vec.at(i)),&(cusum_vec.at(i)),&(ttest_vec.at(i)));
      }
      else{
        kinkidx = FindKinkIdx(&(dirs_vec.at(i)));
      }

      // Don't accept invalid kink indices
      if (kinkidx==-9999) continue;

      // Dont accept kinks too close to the end of the track (let's say within 5 spacepoints)
      if (kinkidx<5 || (dirs_vec.at(i).size()-kinkidx)<5) continue;

      kink_indices->push_back(kinkidx);
      dirs_vec_tmp.push_back(std::vector<double>(dirs_vec.at(i).begin(),dirs_vec.at(i).begin()+kinkidx));
      dirs_vec_tmp.push_back(std::vector<double>(dirs_vec.at(i).begin()+kinkidx+1,dirs_vec.at(i).end()));
      cusum_vec_tmp.push_back(std::vector<double>(cusum_vec.at(i).begin(),cusum_vec.at(i).begin()+kinkidx));
      cusum_vec_tmp.push_back(std::vector<double>(cusum_vec.at(i).begin()+kinkidx+1,cusum_vec.at(i).end()));
      ttest_vec_tmp.push_back(std::vector<double>(ttest_vec.at(i).begin(),ttest_vec.at(i).begin()+kinkidx));
      ttest_vec_tmp.push_back(std::vector<double>(ttest_vec.at(i).begin()+kinkidx+1,ttest_vec.at(i).end()));
    }

    // if kink_indices has the same size as it was before we did the loop above, we found no kinks and should quit
    if (kink_indices->size()==kink_indices_size_before){
      nrpts=10;
      break;
    }

    // If not, set vectors and go for a repeat!
    dirs_vec.clear();
    cusum_vec.clear();
    ttest_vec.clear();
    for (int ivec=0; ivec<dirs_vec_tmp.size(); ivec++){
      dirs_vec.push_back(dirs_vec_tmp.at(ivec));
      cusum_vec.push_back(cusum_vec_tmp.at(ivec));
      ttest_vec.push_back(ttest_vec_tmp.at(ivec));
    }

    nrpts++;
  }

  return kink_indices;

};

#endif
