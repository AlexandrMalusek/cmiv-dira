__kernel void Backproject(const int N, const double ctr, const int numAngles,
                          const int projection_length, __global double *thetaPtr,
                          __global double *projections, __global double *img)
{ 

  int k, y;
  int input_row, output_row;
  
  double xleft = -ctr;
  double ytop  = ctr;
  int center = projection_length/2;

  int x = get_global_id(0);
  double xcoord = xleft + x;

  double cos_theta, sin_theta;
  double t, fraction;
  int a;

  if(x <= N)
  {
    output_row = x*N;

    for(k=0;k<numAngles;++k)
    {
      cos_theta = cos(thetaPtr[k]);
      sin_theta = sin(thetaPtr[k]);

      input_row = k * projection_length;
      t = xcoord * cos_theta + ytop * sin_theta;

      for(y=0;y<N;++y)
      {
        a  = ((int) (t + N)) - N;
        fraction = t - a;
        a +=center;

        img[output_row + y] += fraction * (projections[input_row + a + 1] - projections[input_row + a]) +
                               projections[input_row + a];
        t -= sin_theta;
      }
    }
  }
}
