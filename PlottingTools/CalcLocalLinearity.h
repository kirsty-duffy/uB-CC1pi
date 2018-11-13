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



// Before we can do local linearity we have to make sure the spacepoints are correctly ordered
// A simple approach: start with the spacepoint closest to the track start (or reconstructed neutrino vertex?) and add the closest hit that hasn't been added yet until we have them all
std::vector<std::pair<std::vector<double>,int>> OrderSpacepoints(std::vector<std::vector<double>> spacepoints/*, std::vector<double> startpoint*/){
  std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints;
  // Check start point is valid
  // if (startpoint.at(0)==-999 && startpoint.at(1)==-999 && startpoint.at(2)==-999){
  //   std::cout << "Error: cannot order spacepoints without a valid starting point. Returning an empty vector." << std::endl;
  //   return ordered_spacepoints;
  // }

  // Check input vector contains some points
  if (spacepoints.size()==0){
    std::cout << "Error: spacepoints vector is empty. Returning an empty vector." << std::endl;
    return ordered_spacepoints;
  }

  // Find the closest spacepoint to the start point
  /*double dist = 9999;
  int firstpoint_idx = -9999;
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
  } // end loop to find closest spacepoint to the start*/

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
  int firstpoint_idx = std::distance(dist.begin(), std::max_element(dist.begin(), dist.end()));

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

#endif
