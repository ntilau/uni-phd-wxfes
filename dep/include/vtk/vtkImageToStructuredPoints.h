/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageToStructuredPoints.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkImageToStructuredPoints
 * @brief   Attaches image pipeline to VTK.
 *
 * vtkImageToStructuredPoints changes an image cache format to
 * a structured points dataset.  It takes an Input plus an optional
 * VectorInput. The VectorInput converts the RGB scalar components
 * of the VectorInput to vector pointdata attributes. This filter
 * will try to reference count the data but in some cases it must
 * make a copy.
 */

#ifndef vtkImageToStructuredPoints_h
#define vtkImageToStructuredPoints_h

#include "vtkCommonExecutionModelModule.h" // For export macro
#include "vtkImageAlgorithm.h"

class vtkImageData;
class vtkStructuredPoints;

class VTKCOMMONEXECUTIONMODEL_EXPORT vtkImageToStructuredPoints : public vtkImageAlgorithm
{
public:
  static vtkImageToStructuredPoints* New();
  vtkTypeMacro(vtkImageToStructuredPoints, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Set/Get the input object from the image pipeline.
   */
  void SetVectorInputData(vtkImageData* input);
  vtkImageData* GetVectorInput();
  //@}

  /**
   * Get the output of the filter.
   */
  vtkStructuredPoints* GetStructuredPointsOutput();

protected:
  vtkImageToStructuredPoints();
  ~vtkImageToStructuredPoints() override;

  // to translate the wholeExtent to have min 0 ( I do not like this hack).
  int Translate[3];

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  int FillOutputPortInformation(int, vtkInformation*) override;
  int FillInputPortInformation(int, vtkInformation*) override;

private:
  vtkImageToStructuredPoints(const vtkImageToStructuredPoints&) = delete;
  void operator=(const vtkImageToStructuredPoints&) = delete;
};

#endif
