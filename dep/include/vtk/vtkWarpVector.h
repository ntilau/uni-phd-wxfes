/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkWarpVector.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkWarpVector
 * @brief   deform geometry with vector data
 *
 * vtkWarpVector is a filter that modifies point coordinates by moving
 * points along vector times the scale factor. Useful for showing flow
 * profiles or mechanical deformation.
 *
 * The filter passes both its point data and cell data to its output.
 */

#ifndef vtkWarpVector_h
#define vtkWarpVector_h

#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkPointSetAlgorithm.h"

class VTKFILTERSGENERAL_EXPORT vtkWarpVector : public vtkPointSetAlgorithm
{
public:
  static vtkWarpVector* New();
  vtkTypeMacro(vtkWarpVector, vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Specify value to scale displacement.
   */
  vtkSetMacro(ScaleFactor, double);
  vtkGetMacro(ScaleFactor, double);
  //@}

  int FillInputPortInformation(int port, vtkInformation* info) override;

protected:
  vtkWarpVector();
  ~vtkWarpVector() override;

  int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  double ScaleFactor;

private:
  vtkWarpVector(const vtkWarpVector&) = delete;
  void operator=(const vtkWarpVector&) = delete;
};

#endif
