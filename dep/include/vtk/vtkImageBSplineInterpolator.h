/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageBSplineInterpolator.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkImageBSplineInterpolator
 * @brief   perform b-spline interpolation on images
 *
 * vtkImageBSplineInterpolator can be used to perform b-spline interpolation
 * on images that have been filtered with vtkImageBSplineCoefficients.  The
 * b-spline interpolants provide the maximum possible degree of continuity
 * for a given kernel size, but require that the image data be pre-filtered
 * to generate b-spline coefficients before the interpolation is performed.
 * Interpolating data that has not been pre-filtered will give incorrect
 * results.
 * @sa
 * vtkImageReslice vtkImageBSplineCoefficients vtkBSplineTransform
 * @par Thanks:
 * This class was written by David Gobbi at the Seaman Family MR Research
 * Centre, Foothills Medical Centre, Calgary, Alberta.
 * DG Gobbi and YP Starreveld,
 * "Uniform B-Splines for the VTK Imaging Pipeline,"
 * VTK Journal, 2011,
 * http://hdl.handle.net/10380/3252
 */

#ifndef vtkImageBSplineInterpolator_h
#define vtkImageBSplineInterpolator_h

#include "vtkAbstractImageInterpolator.h"
#include "vtkImagingCoreModule.h" // For export macro

#define VTK_IMAGE_BSPLINE_DEGREE_MAX 9

class vtkImageData;
struct vtkInterpolationInfo;

class VTKIMAGINGCORE_EXPORT vtkImageBSplineInterpolator : public vtkAbstractImageInterpolator
{
public:
  static vtkImageBSplineInterpolator* New();
  vtkTypeMacro(vtkImageBSplineInterpolator, vtkAbstractImageInterpolator);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Set the degree of the spline polynomial.  The default value is 3,
   * and the maximum is 9.  The data must be pre-filtered for the same
   * degree of polynomial with vtkImageBSplineCoefficients.
   */
  void SetSplineDegree(int degree);
  int GetSplineDegree() { return this->SplineDegree; }
  int GetSplineDegreeMinValue() { return 0; }
  int GetSplineDegreeMaxValue() { return VTK_IMAGE_BSPLINE_DEGREE_MAX; }
  //@}

  /**
   * Get the support size for use in computing update extents.  If the data
   * will be sampled on a regular grid, then pass a matrix describing the
   * structured coordinate transformation between the output and the input.
   * Otherwise, pass nullptr as the matrix to retrieve the full kernel size.
   */
  void ComputeSupportSize(const double matrix[16], int support[3]) override;

  /**
   * Returns true if the interpolator supports weight precomputation.
   * This will always return true for this interpolator.
   */
  bool IsSeparable() override;

  //@{
  /**
   * If the data is going to be sampled on a regular grid, then the
   * interpolation weights can be precomputed.  A matrix must be
   * supplied that provides a transformation between the provided
   * extent and the structured coordinates of the input.  This
   * matrix must perform only permutations, scales, and translation,
   * i.e. each of the three columns must have only one non-zero value.
   * A new extent is provided for out-of-bounds checks.
   * THIS METHOD IS THREAD SAFE.
   */
  void PrecomputeWeightsForExtent(const double matrix[16], const int extent[6], int newExtent[6],
    vtkInterpolationWeights*& weights) override;
  void PrecomputeWeightsForExtent(const float matrix[16], const int extent[6], int newExtent[6],
    vtkInterpolationWeights*& weights) override;
  //@}

  /**
   * Free the precomputed weights.  THIS METHOD IS THREAD SAFE.
   */
  void FreePrecomputedWeights(vtkInterpolationWeights*& weights) override;

protected:
  vtkImageBSplineInterpolator();
  ~vtkImageBSplineInterpolator() override;

  /**
   * Update the interpolator.
   */
  void InternalUpdate() override;

  /**
   * Copy the interpolator.
   */
  void InternalDeepCopy(vtkAbstractImageInterpolator* obj) override;

  //@{
  /**
   * Get the interpolation functions.
   */
  void GetInterpolationFunc(
    void (**doublefunc)(vtkInterpolationInfo*, const double[3], double*)) override;
  void GetInterpolationFunc(
    void (**floatfunc)(vtkInterpolationInfo*, const float[3], float*)) override;
  //@}

  //@{
  /**
   * Get the row interpolation functions.
   */
  void GetRowInterpolationFunc(
    void (**doublefunc)(vtkInterpolationWeights*, int, int, int, double*, int)) override;
  void GetRowInterpolationFunc(
    void (**floatfunc)(vtkInterpolationWeights*, int, int, int, float*, int)) override;
  //@}

  /**
   * Build the lookup tables used for the interpolation.
   */
  virtual void BuildKernelLookupTable();

  /**
   * Free the kernel lookup tables.
   */
  virtual void FreeKernelLookupTable();

  int SplineDegree;
  float* KernelLookupTable;

private:
  vtkImageBSplineInterpolator(const vtkImageBSplineInterpolator&) = delete;
  void operator=(const vtkImageBSplineInterpolator&) = delete;
};

#endif
