/**
 * particle_filter.cpp
 *
 * Created on: 20/03/2020
 * Author: Olaoluwa Popoola
 */

//Inititalization section -libraries used in the project
#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using namespace std;

//Initialization of particles
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  num_particles = 500;  //number of particles used for filtering
  normal_distribution<double> dist_x(x, std[0]);  
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  std::default_random_engine gen;

  //particle creation with jitter
  for (int i = 0; i < num_particles; i++) {
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1;    
    particles.push_back(particle);
    weights.push_back(1);
  }
  is_initialized = true; //set initialization to true 
}
    


void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  default_random_engine gen;
  for (int i=0; i<num_particles;i++)
  {
    double new_x, new_y, new_theta;
    //to prevent division by zero
    if(yaw_rate == 0)
    {
      new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
      new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
      new_theta = particles[i].theta;
    }
    //particle prediction using the bicycle model
    else
    {
     new_x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta));
     new_y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
     new_theta = particles[i].theta + yaw_rate*delta_t;
    }
    //adding noise to the prediction 
    normal_distribution<double> dist_x(new_x, std_pos[0]);
    normal_distribution<double> dist_y(new_y, std_pos[1]);
    normal_distribution<double> dist_theta(new_theta, std_pos[2]);
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);                                   
  }

}
                                       

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}
//sc
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */	
  //initializing updates for all particles
  for (int p = 0; p < num_particles; p++) { 
    // initialize sensing variables
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y; 
  
  //get observations
    vector<LandmarkObs> trans_observations;
    LandmarkObs obs;

    //homogeneous transformations using  translation and rotation for the observations
    for(int i=0; i<observations.size(); i++)
    {
      LandmarkObs trans_obs;
      obs = observations[i];
      trans_obs.x = particles[p].x+(obs.x*cos(particles[p].theta)-obs.y*sin(particles[p].theta)); 
      trans_obs.y = particles[p].y+(obs.x*sin(particles[p].theta)+obs.y*cos(particles[p].theta));
      trans_observations.push_back(trans_obs);
    }
    particles[p].weight = 1.0;
    for(int i=0; i<trans_observations.size(); i++) //141
    {
      double closet_dist = sensor_range;
      int association = 0;

      for(int j = 0; j<map_landmarks.landmark_list.size(); j++)//146
      {
        double landmark_x = map_landmarks.landmark_list[j].x_f;
        double landmark_y = map_landmarks.landmark_list[j].y_f;

        double calc_dist = sqrt(pow(trans_observations[i].x - landmark_x,2.0)+pow(trans_observations[i].y - landmark_y,2.0)); 
        if(calc_dist < closet_dist)
        {
          closet_dist = calc_dist;
          association = j;
        }
      }
      if(association != 0)
      {
        double meas_x = trans_observations[i].x;//162
        double meas_y = trans_observations[i].y;
        double mu_x = map_landmarks.landmark_list[association].x_f;
        double mu_y = map_landmarks.landmark_list[association].y_f; 
        //calculate particle weights and multiplier
        long double multiplier = 1/(2*M_PI*std_landmark[0]*std_landmark[1])*exp(-(pow(meas_x-mu_x,2.0)/(2*std_landmark[0]*std_landmark[0])+pow(meas_y-mu_y,2.0)/(2*std_landmark[1]*std_landmark[1])));
        if (multiplier>0)
        {
          particles[p].weight*=multiplier;
        }
      } 
        associations.push_back(association+1);
        sense_x.push_back(trans_observations[i].x);
        sense_y.push_back(trans_observations[i].y);
      
    }
    
       particles[p] = SetAssociations(particles[p], associations, sense_x, sense_y); 
        weights[p] = particles[p].weight; 
  }
}
  
                       


void ParticleFilter::resample() {//#
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
default_random_engine gen;
  discrete_distribution<int> distribution(weights.begin(), weights.end());
  vector<Particle> resample_particles;
  for(int i=0; i<num_particles; i++)
  {
    resample_particles.push_back(particles[distribution(gen)]);
  }
  particles =resample_particles;
}//end chk

Particle ParticleFilter::SetAssociations(Particle particle, 
                                     std::vector<int> associations, 
                                     std::vector<double> sense_x, 
                                     std::vector<double> sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();
  
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
  return particle;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}