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
#include <cfloat>
#include <random> // Need this for sampling from distributions

#include "helper_functions.h"
#include "multiv_gauss.h"

using std::string;
using std::vector;
using std::normal_distribution;

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
   num_particles = 200;  // TODO: Set the number of particles

  // This line creates a normal (Gaussian) distribution for x
   normal_distribution<double> dist_x(x, std[0]);
   normal_distribution<double> dist_y(y, std[1]);
   normal_distribution<double> dist_theta(theta, std[2]);

   for (int i = 0; i < num_particles; ++i) {
	   Particle p;
	   p.id = i;
	   p.x = dist_x(gen);
	   p.y = dist_y(gen);
	   p.theta = dist_theta(gen);
	   p.weight = 1.0;
	   particles.push_back(p);
   }
   is_initialized = true;

   std::cout << "Init Complete!!!" << std::endl;
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
	std::default_random_engine gen;
	double x,y,theta;
	for(int i = 0; i < num_particles; i++)
	{
		if(fabs(yaw_rate) > 0.0001) {
			x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta +(yaw_rate*delta_t))-sin(particles[i].theta));
			y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+(yaw_rate*delta_t)));
			theta = particles[i].theta + (yaw_rate*delta_t);
		}
		else {
			x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
		}

		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_theta(theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
	//std::cout << "Prediction Made!!!" << std::endl;
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
// for each predicated associate it with one of the observations id.
	//std::cout << "Entered dataAssociation!!! p_num:" <<predicted.size() << " o_num:" << observations.size() << std::endl;
	for(unsigned int i = 0; i < observations.size(); i++) {
		int mapId = -1;
		double distance = DBL_MAX;
		for(unsigned int j = 0; j < predicted.size(); j++) {
			double new_distance = dist(observations[i].x,observations[i].y, predicted[j].x, predicted[j].y);
			if(new_distance < distance) {
				distance = new_distance;
				mapId = predicted[j].id;
			}
		}
		observations[i].id = mapId;
	}
	//std::cout << "Exit dataAssociation!!!" << std::endl;
}

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
	// 1. transform observations from Vehicle coordinates to Map coordinates.
	// -------- Homogenous Transformation
	// 2. Associate observations with landmarks using Euclidean Distance.
	// 3. Calulate each particle weight using Multivariate-Gaussian
//for each particle:
//	prob = 1.0
//	for each map_landmarks:
//		dist = particle.x,particle.y,map_landmark.x,landmark.y);
//		prob *= gaussian(dist, sense_noise, observation[i])
//	particle.weight = prob;
	//std::cout << "Entered updateWeights!!! o_num:" << observations.size() << " l_num:" << map_landmarks.landmark_list.size() <<std::endl;

	for(int i = 0; i < num_particles; i++){

		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		vector<LandmarkObs> p;
		for(unsigned int l = 0; l < map_landmarks.landmark_list.size(); l++) {
			double l_x = map_landmarks.landmark_list[l].x_f;
			double l_y = map_landmarks.landmark_list[l].y_f;
			int l_i = map_landmarks.landmark_list[l].id_i;

			double distance = dist(p_x,p_y,l_x,l_y);
			if(distance <= sensor_range) {
				p.push_back( LandmarkObs { l_i, l_x, l_y });
			}
		}

		vector<LandmarkObs> t_obs;
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;

		for(unsigned int o = 0; o < observations.size(); o++) {
			  // transform to map x coordinate
			double x_map = p_x + (cos(p_theta) * observations[o].x) - (sin(p_theta) * observations[o].y);

			// transform to map y coordinate
			double y_map = p_y + (sin(p_theta) * observations[o].x) + (cos(p_theta) * observations[o].y);

			t_obs.push_back( LandmarkObs {observations[o].id, x_map, y_map} );
		}

		dataAssociation(p,t_obs);

		for(unsigned int o = 0; o < t_obs.size(); o++) {
			associations.push_back(t_obs[o].id);
			sense_x.push_back(t_obs[o].x);
			sense_y.push_back(t_obs[o].y);
		}
		SetAssociations(particles[i], associations, sense_x, sense_y);

		particles[i].weight = 1.0;

		for(unsigned int o = 0; o < t_obs.size(); o++) {
			double sig_x, sig_y, x_obs, y_obs, mu_x, mu_y;
			sig_x = std_landmark[0];
			sig_y = std_landmark[1];
			x_obs = t_obs[o].x;
			y_obs = t_obs[o].y;


			for(unsigned int t = 0; t < p.size(); t++) {
				if(p[t].id == t_obs[o].id) {
					mu_x = p[t].x;
					mu_y = p[t].y;
				}
			}

			//std::cout << "sig_x:" << sig_x << " sig_y: " << sig_y << " x_obs: " <<  x_obs << " y_obs: " << y_obs << " mu_x: " <<  mu_x << " mu_y: " << mu_y << std::endl;
			double weight = multiv_prob(sig_x, sig_y, x_obs, y_obs, mu_x, mu_y);

			//std::cout << "weight = " << weight << std::endl;
			particles[i].weight *= weight;
		}

		//std::cout << "particle[ " << i << " ].weight = " << particles[i].weight << std::endl;
	}
	//std::cout << "Exit updateWeights!!!" << std::endl;

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   *//*
	p3 = []
	    index = int(random.random() * N)
	    beta = 0.0
	    mw = max(w)
	    for i in range(N):
	        beta += random.random() * 2.0 * mw
	        while beta > w[index]:
	            beta -= w[index]
	            index = (index + 1) % N
	        p3.append(p[index])*/

	vector<Particle> resampledParticles;
	vector<double> weights;

	//std::cout << "Entered resample!!!" << std::endl;

	for(int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}
	std::default_random_engine gen;
	std::uniform_int_distribution<int> iDist(0, num_particles-1);
	int index = iDist(gen);
	double beta = 0.0;
	int max_index = std::minmax_element(weights.begin(),weights.end()).second - weights.begin();
	double mw = weights[max_index];
	std::uniform_real_distribution<double> rDist(0.0, mw);

	for(int i = 0; i < num_particles; i++) {
		beta += rDist(gen) * 2.0;
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index+1)%num_particles;
		}
		resampledParticles.push_back(particles[index]);
	}

	particles = resampledParticles;

	//std::cout << "Exit resample!!!" << std::endl;

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
