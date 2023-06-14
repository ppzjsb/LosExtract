/********************************************************************************/

/**! \file utils.c 
 *
 * \brief Some general purpose utility functions.
 *
 */

/********************************************************************************/


/** \brief returns the maximum of two doubles            
 *                                                                              
 *   \param x first double
 *   \param y second double
 *
 *   \return maximum of x or y
 */

double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}



/** \brief returns the minimum of two doubles            
 *                                                                              
 *   \param x first double
 *   \param y second double
 *
 *   \return minimum of x or y
 */

double dmin(double x, double y)
{
  if(x < y)
    return x;
  else
    return y;
}


/** \brief linearly interpolates between two doubles
 *                                                                              
 *   \param x  desired abscissa value
 *   \param x0 start of abscissa interval
 *   \param x1 end of abscissa interval
 *   \param y0 start of ordinate interval
 *   \param y1 end of ordinate interval
 *
 *   \return y desired ordinate value
 */

double lerp(double x, double x0, double x1, double y0, double y1)
{
  double y;

  y = y0 + (x-x0) * (y1-y0) / (x1-x0);

  return y;
}
