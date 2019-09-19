/*
 * multiv_gauss.cpp
 *
 *  Created on: Sep 16, 2019
 *      Author: eo144
 */




#include "multiv_gauss.h"
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <iostream>

double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  double exponent;
  double weight;

  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

  // calculate weight using normalization terms and exponent
  weight = gauss_norm * exp(-exponent);

  //std::cout << "norm=" << gauss_norm << " exp=" << exponent  << " w=" << weight << std::endl;
  return weight;
}


#if 0

double homogenous_transform(double x_part,double  y_part,double  x_obs,double  y_obs,double  theta) {
/*
  x_part = 4;
  y_part = 5;
  x_obs = 2;
  y_obs = 2;
  theta = -M_PI/2; // -90 degrees
*/

  // transform to map x coordinate
  double x_map;
  x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);

  // transform to map y coordinate
  double y_map;
  y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);

  return 0;
}

double euclidean_distance(double tobs_x, double tobs_y, double totLandmarks, double** landmarks) {
	double dist = DBL_MIN;

	for(int i = 0; i < totLandmarks; i++)
	{
		dist = std::min(dist, sqrt( pow(tobs_x - landmarks[i][0], 2) + pow(tobs_y - landmarks[i][1], 2) ) );
	}
	return dist;
}
#endif
