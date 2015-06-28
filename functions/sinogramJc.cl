#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

__kernel void sinogramJ(const int no_pixels, __global double *image,
                        __global int *indices, __global int *xdistance,
                        __global int *ydistance, __global double *cosine,
                        __global double *sine, const int xOrigin,
                        __global double *slope, const int M, __global double *projections)
{ 

  int i;
  double pixelvalue;
  double r, fraction;
  int r_index;
  double distance;
  double leftdistance, rightdistance;
  double slopedpixelvalue, leftpixel, rightpixel;

  double max_r = 0;

  int k = get_global_id(0);

  for(i=0;i<no_pixels;++i)
  {
    pixelvalue = image[indices[i]];

    /* Find the index for the radial coordinates */
    r = xdistance[i]*cosine[k] + ydistance[i]*sine[k];      

    /* add xOrigin to shift center of image, avoiding negative values */
    r += xOrigin;   
    r_index = (int) r;  
    fraction = r - r_index;

    /* Get the pixel value and distribute between two pixels
     * Calculates the distance once as it is used multiple times
     * The slope is dependent on the angle, decreasing the
     * triangle size */
    distance = fraction*slope[k];

    /* No contribution if the distance is less than 0 */
    /* Equal to 
     * (1 - fraction*slope[k]) and
     * (1 - (1 - fraction) * slope[k])*/
    leftdistance  = MAX(0, (1 - distance));
    rightdistance = MAX(0, (1 + distance - slope[k]));

    slopedpixelvalue = pixelvalue * slope[k];
    leftpixel  = leftdistance  * slopedpixelvalue;
    rightpixel = rightdistance * slopedpixelvalue;

    projections[k*M + r_index + 0] += leftpixel;
    projections[k*M + r_index + 1] += rightpixel;

  }
}
