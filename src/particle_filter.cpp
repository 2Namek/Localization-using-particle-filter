/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

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

double weights_sum=0.0;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
 /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
   std::default_random_engine gen;
 normal_distribution<double> dist_x(x, std[0]);
 normal_distribution<double> dist_y(y, std[1]);
 normal_distribution<double> dist_theta(theta, std[2]);
   
  num_particles = 20;// TODO: Set the number of particles
  for(int i=0;i<num_particles;i++)
  {
		  // double sample_x, sample_y, sample_theta;
		  // sample_x=dist_x(gen);
		  // sample_y=dist_y(gen);
		  // sample_theta=dist_theta(gen);
		  // particlesx.push_back(sample_x);
		  // particlesy.push_back(sample_y);
		  // particlestheta.push_back(sample_theta);
		  // weights.push_back(1);
	
		  Particle p;
		  p.x=dist_x(gen);
		  p.y=dist_y(gen);
		  p.theta=dist_theta(gen);
		  p.weight=1;
		  particles.push_back(p);
		  weights.push_back(1);
		  
		  
  }
  
 //uncomment to see the output of this function
  for(int i=0;i<num_particles;i++)
  {
	std::cout<<"initparticlex="<<particles[i].x<<std::endl;
	std::cout<<"initparticley="<<particles[i].y<<std::endl;
	std::cout<<"initparticletheta="<<particles[i].theta<<std::endl;
	std::cout<<"initweights"<<particles[i].weight<<std::endl;
  }
  
  is_initialized = true;
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
   // std::cout<<"start pred"<<std::endl;
    std::default_random_engine gen;
 	  normal_distribution<double> dist_x(0, std_pos[0]);
      normal_distribution<double> dist_y(0, std_pos[1]);
      normal_distribution<double> dist_theta(0, std_pos[2]);
   for(int i=0;i<num_particles;i++)
  {
	if (fabs(yaw_rate) > 0.00001) {
	  particles[i].x=particles[i].x+((velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta)))+dist_x(gen);
	  particles[i].y=particles[i].y+((velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t)))+dist_y(gen);
	  particles[i].theta=particles[i].theta+(yaw_rate*delta_t)+dist_theta(gen);
	  } else {
		  particles[i].x  = particles[i].x  + velocity * delta_t * cos(particles[i].theta)+dist_x(gen);
			particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta)+dist_y(gen);
			 particles[i].theta= particles[i].theta+dist_theta(gen);
	  }
  }
    // for(int i=0;i<num_particles;i++)
  // {
	// std::cout<<"predparticlex="<<particles[i].x<<std::endl;
	// std::cout<<"predparticley="<<particles[i].y<<std::endl;
	// std::cout<<"predparticletheta="<<particles[i].theta<<std::endl;
	// std::cout<<"predweights"<<particles[i].weight<<std::endl;
  // }
     // std::cout<<"end pred"<<std::endl;
 

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

double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
//std::cout<<"gauss_norm= "<<gauss_norm<<std::endl;
  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
  //std::cout<<"exponent= "<<exponent<<std::endl;  
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
  //std::cout<<"weight= "<<weight<<std::endl;    
  return weight;
}
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
    /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can rea
   d more about this distribution here: 
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
   double tempx;
   double tempy;
   weights_sum = 0.0;
	double highest_weight = -1.0;
	int best_particle;
	for(int i=0;i<num_particles;i++)
  {
		// std::cout<<"particlex="<<particles[i].x<<std::endl;
		// std::cout<<"particley="<<particles[i].y<<std::endl;
		// std::cout<<"particletheta="<<particles[i].theta<<std::endl;
		//particles[i].weight = 1.0;
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();
		
	  for(unsigned int j=0;j<observations.size();j++)
	  {
		   // std::cout<<"observations[j].x= "<<observations[j].x<<std::endl;
		   // if (((observations[j].x<sensor_range) && (observations[j].x>-sensor_range)) && ((observations[j].y<sensor_range) && (observations[j].y>-sensor_range)))
		   // {
	tempx=particles[i].x+ (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
	particles[i].sense_x.push_back(tempx);
		// std::cout<<"observations[j].y= "<<observations[j].y<<std::endl;
	tempy=particles[i].y+ (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
	particles[i].sense_y.push_back(tempy);

		// std::cout<<"particle_sense_x="<<particles[i].sense_x[j]<<"j="<<j<<std::endl;
		// std::cout<<"particle_sense_y="<<particles[i].sense_y[j]<<"j="<<j<<std::endl;		
	  // }
	  }
	  
	
	    particles[i].associations.clear();
		particles[i].associationss.clear();
	  for(unsigned int q=0;q<particles[i].sense_x.size();++q)
  {
	int closest_landmark = 0;
    int min_dist = 999999;
    int curr_dist;
    // Iterate through all landmarks to check which is closest
	  
	
    for (unsigned int z = 0; z < map_landmarks.landmark_list.size(); ++z) 
	{
      // Calculate Euclidean distance
       curr_dist = sqrt(pow(particles[i].sense_x[q] - map_landmarks.landmark_list[z].x_f, 2)
                     + pow(particles[i].sense_y[q] - map_landmarks.landmark_list[z].y_f, 2));
	 // curr_dist=dist(map_landmarks.landmark_list[z].x_f,map_landmarks.landmark_list[z].y_f,particles[i].sense_x[q],particles[i].sense_y[q]);
	  // std::cout<<"z="<<z<<" q="<<q<<" curr_dist="<<curr_dist<<std::endl;
      // Compare to min_dist and update if closest
      if (curr_dist < min_dist) 
	  { 
        min_dist = curr_dist;
        closest_landmark = z;
		// std::cout<<"min_dist="<<min_dist<<"closest_landmark="<<closest_landmark<<std::endl; 
		
      }
    }
	if(min_dist<sensor_range)
	{
	 particles[i].associations.push_back(map_landmarks.landmark_list[closest_landmark].id_i); 
	 particles[i].associationss.push_back(closest_landmark);
	 std::cout<<"id_i="<<map_landmarks.landmark_list[closest_landmark].id_i<<"closest_landmark="<<closest_landmark<<std::endl;
  }
  else
  {
	  
       particles[i].sense_x.erase(particles[i].sense_x.begin() + q);
	  particles[i].sense_y.erase(particles[i].sense_y.begin() + q);
  }
  }
	  // for (unsigned int q=0;q<particles[i].sense_x.size();++q)
  // {
	
	// std::cout<<"associations="<<particles[i].associations[q]<<std::endl;
  // }
	  
	  
	  
	  
	  
	  
  }
  
  // weights.clear();
  	for(int i=0;i<num_particles;i++)
 {
	 // int previousweight=1;
	 // double tempweight[observations.size()];
	 double multiplication=1;
  for (unsigned int q=0;q<particles[i].sense_x.size();++q)
  {
	  int tempid=particles[i].associationss[q];
	  // std::cout<<"tempid= "<<tempid<<std::endl;
	  // std::cout<<"std_landmark[0]= "<<std_landmark[0]<<std::endl;
	  // std::cout<<"std_landmark[1]= "<<std_landmark[1]<<std::endl;
	  // std::cout<<"particles[i].sense_x[q]= "<<particles[i].sense_x[q]<<std::endl;
	  // std::cout<<"particles[i].sense_y[q]= "<<particles[i].sense_y[q]<<std::endl;
	  // std::cout<<"map_landmarks.landmark_list[tempid].x_f= "<<map_landmarks.landmark_list[tempid].x_f<<std::endl;
	  // std::cout<<"map_landmarks.landmark_list[tempid].y_f= "<<map_landmarks.landmark_list[tempid].y_f<<std::endl;
	  double tempweight = multiv_prob(std_landmark[0],std_landmark[1],particles[i].sense_x[q],particles[i].sense_y[q],map_landmarks.landmark_list[tempid].x_f,map_landmarks.landmark_list[tempid].y_f);
	  // std::cout<<"tempweight= "<<tempweight<<"  "<<q<<std::endl;
	   // if (tempweight > 0.00001) 
	   // {
	  // particles[i].weight=particles[i].weight*tempweight;
	  multiplication=multiplication*tempweight;
	  
	   // }
	   // else
	   // {
		  // multiplication=0.00000001; 
	   // }
	  
	  
  }
  particles[i].weight=multiplication;
weights[i]=particles[i].weight;
  //std::cout<<"multiplication= "<<multiplication<<std::endl;
  
 
 }

 // std::cout<<"test0"<<std::endl;

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
	std::vector<Particle> resampled_particles;

	double max_weight = *std::max_element(weights.begin(),weights.end());

	std::default_random_engine gen;
	std::uniform_real_distribution<double> random_distr(0.0, max_weight);
	std::uniform_int_distribution<int> index_random_distr(0.0, (num_particles-1));

	int index = index_random_distr(gen);
	double beta = 0.0;


	resampled_particles.clear();

	for (int i = 0; i < num_particles; ++i)
	{
		beta += 2.0 * random_distr(gen);
		while (beta > (weights[index]/*/weights_sum*/))
		{
			beta -= (weights[index]/*/weights_sum*/);
			index = (index + 1) % num_particles;
		}

		resampled_particles.push_back(particles[index]);
	}
   particles = resampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
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
