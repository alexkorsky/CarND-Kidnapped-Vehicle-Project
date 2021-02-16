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

static std::default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 500;  // TODO: Set the number of particles


  // Normal distributions
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  // Generate particles with normal distribution with mean on GPS values.
  for (int i = 0; i < num_particles; i++)
  {
      Particle p;
      p.id = i;
      p.x = dist_x(gen);
      p.y = dist_y(gen);
      p.theta = dist_theta(gen);
      p.weight = 1.0;

      particles.push_back(p);
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

	// For each particle we predict its position after delta_t time given the particle's
	// current position and velocity vector and some noise

    //Normal distributions for sensor noise
    std::normal_distribution<double> dist_x(0, std_pos[0]);
    std::normal_distribution<double> dist_y(0, std_pos[1]);
    std::normal_distribution<double> dist_theta(0, std_pos[2]);

    for (int i = 0; i < num_particles; i++)
    {
    	// protecting against very small yaw_rate -- otherwise particles will not be "moving"
        if (fabs(yaw_rate) >= 0.00001)
        {
            particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
            particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
            particles[i].theta += yaw_rate * delta_t;
        }
        else
        {
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
        }

        // Add some noise
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
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

    // O(ij);
    for (unsigned int i = 0; i < observations.size(); i++)
    {
        double min_dist = std::numeric_limits<double>::max();

        int closest_predicted_id = -1;
        for (unsigned int j = 0; j < predicted.size(); j++)
        {

        	// dist() method is defined in helper_funnctions.h
            double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

            if (distance < min_dist)
            {
                min_dist = distance;
                closest_predicted_id = predicted[j].id;
            }
        }

        observations[i].id = closest_predicted_id;
    }

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

	// here we compute the "weight" of each particle -- meaning how confident we are
	// that a particle's location is actually our car location

	// How it is done: simple
	//   -- take observation to landmarks we actually have (with some noise)
	//   -- translate them to what they would be if we translated from our position to particle's position
	//   -- compare the difference of translated observations and particle's own observations
	//   -- calculate the weight based on this difference


    // loop through each particle
    for (int i = 0; i < num_particles; i++)
    {
        double p_x = particles[i].x;
        double p_y = particles[i].y;
        double p_theta = particles[i].theta;

        // Create a vector of Map landmarks that are within the sensor range of this particle

        vector<LandmarkObs> landmarksWithinRange;
        // Each map landmark for loop
        for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
        {
            float landmark_x = map_landmarks.landmark_list[j].x_f;
            float landmark_y = map_landmarks.landmark_list[j].y_f;
            int id = map_landmarks.landmark_list[j].id_i;

            if (fabs(p_x - landmark_x) <= sensor_range && fabs(p_y - landmark_y) <= sensor_range)
            {
            	landmarksWithinRange.push_back(LandmarkObs{id, landmark_x, landmark_y});
            }
        }


        // from each particle's view transform  observations  from vehicle coordinates to map coordinates
        vector<LandmarkObs> transformedObs;

        for (unsigned int j = 0; j < observations.size(); j++)
        {
            double transformed_x = p_x + cos(p_theta) * observations[j].x - sin(p_theta) * observations[j].y;
            double transformed_y = p_y + sin(p_theta) * observations[j].x + cos(p_theta) * observations[j].y;
            transformedObs.push_back(LandmarkObs{observations[j].id, transformed_x, transformed_y});
        }


        // Data association for the predictions and transformed observations on current particle
        dataAssociation(landmarksWithinRange, transformedObs);


        particles[i].weight = 1.0;

        for (unsigned int j = 0; j < transformedObs.size(); j++)
        {
            double transformedObs_x = transformedObs[j].x;
            double transformedObs_y = transformedObs[j].y;
            int associated_landmarkId = transformedObs[j].id;


            // get the x,y coordinates of the prediction associated with the current observation
            double landmark_x, landmark_y;
            for (unsigned int k = 0; k < landmarksWithinRange.size(); k++)
            {
              if (landmarksWithinRange[k].id == associated_landmarkId)
              {
            	  landmark_x = landmarksWithinRange[k].x;
            	  landmark_y = landmarksWithinRange[k].y;

            	  break;
              }
            }

            // Weight for this observation with multivariate Gaussian
            double dX = transformedObs_x - landmark_x;
            double dY = transformedObs_y - landmark_y;

            double std_x = std_landmark[0];
            double std_y = std_landmark[1];

            double weight = (1 / (2 * M_PI * std_x * std_y)) * exp(-(dX * dX / (2 * std_x * std_x) + (dY * dY / (2 * std_y * std_y))));

            particles[i].weight = particles[i].weight * weight;

            /*
            if (weight == 0)
            {
                particles[i].weight = particles[i].weight * 0.00001;

            }
            else
            {
                particles[i].weight = particles[i].weight * weight;
            }
            */
        }
    }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

	  vector<Particle> resampled_particles;

	  // get all of the current weights
	  vector<double> weights;
	  for (int i = 0; i < num_particles; i++)
	  {
		  weights.push_back(particles[i].weight);
	  }

	  // generate random starting index for resampling wheel
	  std::uniform_int_distribution<int> uniint_dist(0, num_particles-1);
	  auto index = uniint_dist(gen);

	  // get max weight
	  double max_weight = *max_element(weights.begin(), weights.end());

	  // uniform random distribution [0.0, max_weight)
	  std::uniform_real_distribution<double> unireal_dist(0.0, max_weight);

	  double beta = 0.0;

	  // spin the resample wheel!
	  for (int i = 0; i < num_particles; i++)
	  {
		  beta += unireal_dist(gen) * 2.0;
		  while (beta > weights[index])
		  {
			  beta -= weights[index];
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
